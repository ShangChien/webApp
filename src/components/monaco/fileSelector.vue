<!-- eslint-disable vue/attributes-order -->
<script setup lang='ts'>
import { computed, inject, onMounted, reactive, ref } from 'vue'
import type { VNodeChild } from 'vue'
import { useDropZone, useElementSize, useFileDialog } from '@vueuse/core'
import { NInput, NPopover, NScrollbar, useMessage } from 'naive-ui'
import editor from '@/components/monaco/editor.vue'
import filehandle from '@/components/monaco/filehandle.vue'

const headerTabName = inject('headerTabName', 'files')
const message = useMessage()
function info(type: any | string, content: string | (() => VNodeChild)) {
  message.create(
    content,
    {
      type,
      closable: true,
      duration: 3000,
      keepAliveOnHover: true,
    },
  )
}
const textFiles = ref<{ name: string; contents: string }[]>([])
const hasFile = computed(() => textFiles.value.length > 0)
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
const { width: _width, height: _height } = useElementSize(dropZoneRef)
const sizewidth = computed(() => `${_width.value + 12}px`)
const sizeheight = computed(() => `${_height.value + 12}px`)
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

const showPop = ref(false)
const fileName = ref('')
const vaildName = computed(() => (!textFiles.value.some(i => i.name === fileName.value)) && (fileName.value !== ''))
function newFile(name: string) {
  if (vaildName.value) {
    textFiles.value.push({ name, contents: '' })
    edit(textFiles.value.length - 1)
  } else {
    info('warning', 'invalid name or duplicated name in files')
  }
  showPop.value = false
  console.log(textFiles.value)
}

function edit(index: number) {
  textEditor.strText = textFiles.value[index].contents
  textEditor.name = textFiles.value[index].name
  textEditor.id = index /* last to set id */
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
  textEditor.id = null
}

onMounted(() => {
})
</script>

<template>
  <div class="h-full w-full flex flex-nowrap items-center bg-slate-1 rd-2 max-w-full py-1 box-border">
    <div
      ref="dropZoneRef" class="flex-none flex-(~ col nowrap) justify-start items-center rd-1 m-1 h-full box-border w-280px bg-green-50"
      :class="[isOverDropZone ? 'outline outline-0.2rem outline-green-4 bg-green-1' : '']"
    >
      <div class="flex flex-nowrap flex-none box-border justify-center items-center mx-1 p-0">
        <NPopover trigger="manual" overlap placement="top" style="padding: 0" :show="showPop" @clickoutside="showPop = false">
          <template #trigger>
            <div
              class="bg-sky-2 rd-1 m-1 mb-0 p-1 text-(l blue-5) leading-1em cursor-pointer
              hover:(bg-sky-3)
              active:(outline outline-2px outline-blue-3)"
              @click="showPop = !showPop"
            >
              New
            </div>
          </template>
          <div class="flex flex-nowrap justify-center items-center">
            <NInput v-model:value="fileName" size="small" placeholder="name" class="m-1" />
            <div
              class="bg-green-2 rd-1 mr-1 p-1 text-(lg green-5) leading-1em cursor-pointer text-nowrap box-border
              hover:(bg-green-3)
              active:(outline outline-2px outline-green-3)"
              @click="newFile(fileName)"
            >
              ok
            </div>
          </div>
        </NPopover>

        <div
          class="bg-sky-2 rd-1 m-1 mb-0 p-1 text-(l blue-5) leading-1em cursor-pointer
          hover:(bg-sky-3)
          active:(outline outline-2px outline-blue-3)"
          @click="open()"
        >
          Open
        </div>
        <div
          v-show="!!files"
          class="bg-sky-2 rd-1 m-1 mb-0 p-1 text-(l blue-5) leading-1em cursor-pointer
          hover:(bg-sky-3)
          active:(outline outline-2px outline-blue-3)"
          @click="clear()"
        >
          Clear
        </div>
      </div>
      <div class="flex-auto maxH w-full box-border p-1 relative bg-green-50">
        <NScrollbar class="max-h-full flex-auto rd-1 bg-lightblue-1 relative">
          <div v-if="hasFile" class="flex-(~ col nowrap) justify-start items-center box-border h-full w-full gap-1 p-1 rd-1 ">
            <div
              v-for="(file, index) in textFiles"
              :key="index"
              class="rd-2 w-full "
              :class="[index === textEditor.id ? 'bg-slate-50' : 'bg-slate-2']"
            >
              <filehandle
                :fileName="file.name"
                :index="index"
                @view="(i) => edit(i)"
                @edit="(i) => edit(i)"
                @delete="(i) => deleteItem(i)"
              />
            </div>
          </div>
          <div v-else class="text-center h-full w-full box-border">
            Click open button to pick files<br>
            Or drag files to here<br>
            (not support directory)
          </div>
        </NScrollbar>
      </div>
    </div>
    <div class="h-full size relative box-border bg-white rd-1">
      <editor
        v-if="headerTabName === 'files'"
        :strText="textEditor.strText"
        :name="textEditor.name"
        class="h-full w-full box-border"
        @sync="(e) => updateText(e)"
      />
      <div
        v-else-if="headerTabName === 'ml'"
      >
        ssss
      </div>
    </div>
  </div>
</template>

<style>
.size {
  width:calc(100% - v-bind(sizewidth))
}
.maxH {
  max-height:calc(v-bind(sizeheight) - 38px)
}
</style>
