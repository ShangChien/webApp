<!-- eslint-disable vue/attributes-order -->
<script setup lang='ts'>
import { h, onMounted, reactive, ref, watch } from 'vue'
import type { VNodeChild } from 'vue'
import { useMessage } from 'naive-ui'
import editor from '@/components/monaco/editor.vue'
import mlPredict from '@/components/ml/mlPredict.vue'
import viewSetting from '@/components/visualAnalysis/viewSetting.vue'
import fileSelector from '@/components/monaco/fileSelector.vue'
import { usemlStore } from '@/components/ml/mlStore'

const fileIndex = ref<number | null>(null)
const allFiles = ref<{ name: string; contents: string }[]>([])

watch(fileIndex, (v) => { console.log(allFiles.value[v]?.name ?? '') })

const MLStore = usemlStore()

const Tabs = reactive({
  files: () => h(editor, {
    strText: allFiles.value[fileIndex.value]?.contents ?? '',
    name: allFiles.value[fileIndex.value]?.name ?? '',
    onSync: e => updateText(e),
  }, null),
  mlPredict: () => h(mlPredict, {
    allFiles: allFiles.value,
    index: fileIndex.value,
  }),
  visualAnalysis: () => h(viewSetting, {
    'allFiles': allFiles.value,
    'index': fileIndex.value,
    'onUpdate:index': (index: number | null) => fileIndex.value = index,
  }),
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
        <component :is="Tabs[MLStore.currentTab]" />
      </KeepAlive>
    </div>
  </div>
</template>

<style>
</style>
