<!-- eslint-disable vue/attributes-order -->
<script setup lang='ts'>
import { computed, h, inject, onMounted, ref, watch } from 'vue'
import type { VNodeChild } from 'vue'
import { useElementSize } from '@vueuse/core'
import { useMessage } from 'naive-ui'
import editor from '@/components/monaco/editor.vue'
import mlFlow from '@/components/ml/mlFlow.vue'
import fileSelector from '@/components/monaco/fileSelector.vue'

const fileIndex = ref<number | null>(null)
const allFiles = ref<{ name: string; contents: string }[]>([])
const _hasFile = computed(() => allFiles.value.length > 0)

watch(fileIndex, (v) => { console.log(allFiles.value[v]?.name ?? '') })

const dom = ref()
const { width: _width } = useElementSize(dom)
const sizewidth = computed(() => `${_width.value + 12}px`)

const headerTabName = inject('headerTabName', 'files') // files or ml
const Tabs = {
  files: () => h(editor, {
    strText: allFiles.value[fileIndex.value]?.contents ?? '',
    name: allFiles.value[fileIndex.value]?.name ?? '',
    onSync: e => updateText(e),
  }, null),
  ml: mlFlow,
}

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
  <div class="h-full w-full flex flex-nowrap justify-start items-center bg-slate-1 rd-2 max-w-full p-1 box-border gap-1">
    <fileSelector
      ref="dom"
      v-model:allFiles="allFiles"
      v-model:fileIndex="fileIndex"
    />
    <div class="flex-auto h-full sizeW relative box-border bg-white rd-1 relative">
      <KeepAlive>
        <component :is="Tabs[headerTabName]" />
      </KeepAlive>
    </div>
  </div>
</template>

<style>
.sizeW {
  width:calc(100% - v-bind(sizewidth))
}
</style>
