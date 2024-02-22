<!-- eslint-disable vue/attributes-order -->
<script setup lang='ts'>
import { h, onMounted, reactive, ref, watch } from 'vue'
import type { VNodeChild } from 'vue'
import { useMessage } from 'naive-ui'
import editor from '@/components/monaco/editor.vue'
import spectrum from '@/components/dataCheck/spectrum.vue'
import fileSelector from '@/components/monaco/fileSelector.vue'
import { useDataCheckStore } from '@/components/dataCheck/dataCheckStore'

const fileIndex = ref<number | null>(null)
const allFiles = ref<{ name: string; contents: string }[]>([])

watch(fileIndex, (v) => { console.log(allFiles.value[v]?.name ?? '') })

const Store = useDataCheckStore()

const Tabs = reactive({
  files: () => h(editor, {
    strText: allFiles.value[fileIndex.value]?.contents ?? '',
    name: allFiles.value[fileIndex.value]?.name ?? '',
    onSync: e => updateText(e),
  }, null),
  spectrum: () => h(spectrum, {
    allFiles: allFiles.value,
    index: fileIndex.value,
  }),
  table: () => h('div', null, ['to do ...']),
})

const message = useMessage()
function _info(type: any | string, content: string | (() => VNodeChild)) {
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

function updateText(newStrText: string) {
  allFiles.value[fileIndex.value].contents = newStrText
}
onMounted(() => {
})
</script>

<template>
  <div class="h-full w-full flex-auto flex flex-nowrap justify-start items-center bg-slate-1 rd-2 max-w-full p-1 box-border gap-1">
    <fileSelector
      v-model:allFiles="allFiles"
      v-model:fileIndex="fileIndex"
    />
    <div class="flex-auto h-full max-h-full relative box-border bg-white rd-1 relative flex">
      <KeepAlive>
        <component :is="Tabs[Store.currentTab]" v-if="fileIndex !== null" />
      </KeepAlive>
    </div>
  </div>
</template>

<style>
</style>
