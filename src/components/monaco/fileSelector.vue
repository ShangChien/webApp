<script setup lang='ts'>
import { computed, onMounted, reactive, ref } from 'vue'
import { useDropZone, useElementSize, useFileDialog } from '@vueuse/core'
import editor from '@/components/monaco/editor.vue'

const textFiles = ref<{ name: string; contents: string }[]>([])
const textEditor = reactive<{ id: number | null; name: string; strText: string }>({
  id: null,
  name: '',
  strText: 'sss',
})
function updateText(newStrText: string) {
  textEditor.strText = newStrText
  textFiles.value[textEditor.id].contents = newStrText
}

const dropZoneRef = ref(null)
const { width: _width } = useElementSize(dropZoneRef)
const sizewidth = computed(() => `${_width.value + 12}px`)
const { isOverDropZone } = useDropZone(dropZoneRef, onDrop)
function onDrop(files: File[]) {
  readfiles(files)
}

const { files, open, reset, onChange } = useFileDialog()
onChange((files) => {
  readfiles(files)
})

function _readfile(file: any) {
  const getdText = new Promise((resolve, reject) => {
    const reader = new FileReader()
    reader.onload = (e) => {
      const contents = e.target.result
      resolve({ name: file.name, contents })
    }
    reader.onerror = function (e) {
      reject(e)
    }
    reader.readAsText(file)
  })
  Promise.resolve(getdText)
    .then((file: { name: string; contents: string }) => {
      textEditor.strText = file.contents
      textEditor.name = file.name
      console.log(`读取${file.name}文件`)
    })
    .catch(e => console.log(`error:${e}`))
}

function readfiles(files: any) {
  const fileContents = Array.from(files).map((file: File) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (e) => {
        const contents = e.target.result
        resolve({ name: file.name, contents })
      }
      reader.onerror = function (e) {
        reject(e)
      }
      reader.readAsText(file)
    })
  })
  Promise.all(fileContents)
    .then((files: { name: string; contents: string }[]) => {
      textFiles.value = files
      console.log(`读取${files?.length}个文件`)
    })
    .catch(e => console.log(`error:${e}`))
}

function deleteItem(index: number) {
  textFiles.value.splice(index, 1)
  if (index === textEditor.id) {
    textEditor.strText = ''
    textEditor.id = null
  }
}

function clear() {
  reset()
  textFiles.value = []
  textEditor.name = ''
  textEditor.strText = ''
  textEditor.name = ''
}

onMounted(() => {
})
</script>

<template>
  <div class="h-full w-full flex flex-nowrap items-center bg-slate-1 rd-2 min-w-280px max-w-full py-1 box-border">
    <div
      ref="dropZoneRef" class="flex-none flex-(~ col nowrap) justify-center items-center rd-1 m-1 h-full bg-green-50"
      :class="[isOverDropZone ? 'outline outline-0.2rem outline-green-4 bg-green-1' : '']"
    >
      <div class="flex flex-nowrap flex-none justify-center items-center m-1 p-0">
        <div
          class="bg-sky-2 rd-1 m-1 p-1 text-(l blue-5) leading-1em cursor-pointer
          hover:(bg-sky-3)
          active:(outline outline-2px outline-blue-3)"
          @click="open()"
        >
          Open
        </div>
        <div
          v-show="!!files"
          class="bg-sky-2 rd-1 m-1 p-1 text-(l blue-5) leading-1em cursor-pointer
          hover:(bg-sky-3)
          active:(outline outline-2px outline-blue-3)"
          @click="clear()"
        >
          Clear
        </div>
      </div>
      <div v-if="files" class="flex-(~ col nowrap) justify-center items-center gap-1 m-1 p-1 rd-2 bg-blue-1">
        <p>You have selected: <b>{{ `${files.length} ${files.length === 1 ? 'file' : 'files'}` }}</b></p>
        <div
          v-for="(file, index) in textFiles"
          :key="file.name"
          :class="[index === textEditor.id ? 'bg-emerald-1' : 'bg-light']"
          class="flex-(~ nowrap none) justify-center items-center rd-2"
        >
          <div class="px-1">
            {{ file.name }}
          </div>
          <div
            class="bg-sky-2 rd-1 m-1 p-1 text-(l blue-5) leading-1em cursor-pointer
            hover:(bg-sky-3)
            active:(outline outline-2px outline-blue-3)"
            @click="() => { textEditor.strText = file.contents; textEditor.name = file.name; textEditor.id = index /* last to set id */ }"
          >
            view
          </div>
          <div
            class="bg-red-2 rd-1 m-1 p-1 text-(l blue-5) leading-1em cursor-pointer
            hover:(bg-red-3)
            active:(outline outline-2px outline-red-3)"
            @click=" deleteItem(index)"
          >
            del
          </div>
        </div>
      </div>
      <div v-else class="text-center m-1">
        Click open button to pick files<br>
        Or drag files to here (not support directory)
      </div>
    </div>
    <div class="h-full size relative box-border">
      <editor
        :strText="textEditor.strText"
        :name="textEditor.name"
        class="h-full w-full box-border"
        @update:strText="(e) => updateText(e)"
      />
    </div>
  </div>
</template>

<style>
.size {
  width:calc(100% - v-bind(sizewidth))
}
</style>
