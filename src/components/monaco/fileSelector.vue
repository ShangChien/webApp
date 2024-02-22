<!-- eslint-disable vue/attributes-order -->
<script setup lang='ts'>
import { computed, ref } from 'vue'
import type { VNodeChild } from 'vue'
import { useDropZone, useElementSize, useFileDialog } from '@vueuse/core'
import { NInput, NPopover, NScrollbar, useMessage } from 'naive-ui'
import filehandle from '@/components/monaco/filehandle.vue'

const fileIndex = defineModel<number | null>('fileIndex')
const allFiles = defineModel<{ name: string; contents: string }[]>('allFiles')
const hasFile = computed(() => allFiles.value.length > 0)

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

const dropZoneRef = ref(null)
const { width: _width, height: _height } = useElementSize(dropZoneRef)
const sizeheight = computed(() => `${_height.value - 32}px`)
const { isOverDropZone } = useDropZone(dropZoneRef, onDrop)
function onDrop(files: File[]) {
  readfiles(files)
}

const { files, open, reset, onChange } = useFileDialog()
onChange((files) => {
  readfiles(files)
})

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
      allFiles.value = files
      fileIndex.value = null
      console.log(`读取${files?.length}个文件`)
    })
    .catch(e => console.log(`error:${e}`))
}

const showPop = ref(false)
const fileName = ref('')
const vaildName = computed(() => (!allFiles.value.some(i => i.name === fileName.value)) && (fileName.value !== ''))
function newFile(name: string) {
  if (vaildName.value) {
    allFiles.value.push({ name, contents: '' })
    edit(allFiles.value.length - 1)
  } else {
    info('warning', 'invalid name or duplicated name in files')
  }
  showPop.value = false
  console.log(allFiles.value)
}

function edit(index: number) {
  fileIndex.value = index
}

function deleteItem(index: number) {
  if (fileIndex.value === index) {
    fileIndex.value = null
  }
  allFiles.value.splice(index, 1)
}

function clear() {
  reset()
  allFiles.value = []
}
</script>

<template>
  <div
    ref="dropZoneRef" class="flex-none flex-(~ col nowrap) justify-start items-center rd-1 h-full box-border max-h-full w-280px bg-green-50"
    :class="[isOverDropZone ? 'outline outline-0.2rem outline-green-4 bg-green-1' : '']"
  >
    <div class="flex flex-nowrap flex-none box-border justify-center items-center p-0">
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
    <div class="flex-auto w-full box-border bg-green-50 relative p-1 rd-1">
      <NScrollbar :style="{ maxHeight: sizeheight }" class="rd-1 w-full box-border bg-sky-1">
        <div v-if="hasFile" class="flex-(~ col nowrap) justify-start items-center box-border w-full gap-1 p-1 rd-1 ">
          <div
            v-for="(file, index) in allFiles"
            :key="index"
            class="rd-2 w-full "
            :class="[index === fileIndex ? 'bg-slate-50' : 'bg-slate-2']"
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
        <div v-else class="text-center w-full box-border">
          Click open button to pick files<br>
          Or drag files to here<br>
          (not support directory)
        </div>
      </NScrollbar>
    </div>
  </div>
</template>

<style>
</style>
