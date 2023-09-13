<script setup lang='ts'>
import { ref } from 'vue'
import MonacoEditor from 'monaco-editor-vue3'
import { useElementSize, watchDebounced } from '@vueuse/core'

const props = defineProps<{ name: string }>()
const strText = defineModel<string>('strText')
const editorDom = ref(null)
const renderKey = ref<number>(0)
const { width: editorW, height: editorH } = useElementSize(editorDom)
const options = {
  colorDecorators: true,
  lineHeight: 2,
  tabSize: 2,
  theme: 'vs',
  language: 'python',
}

watchDebounced(
  [editorH, editorW],
  () => {
    renderKey.value += 1
    console.log('key', editorW.value, editorH.value, renderKey.value)
  },
  { debounce: 1000, maxWait: 5000 },
)
</script>

<template>
  <div ref="editorDom" class="flex-(~ col nowrap) justify-center items-center box-border p-1 bg-sky-50 rd-2">
    <div class="flex-none">
      {{ `file name: ${props.name}` }}
    </div>
    <div :key="renderKey" class="flex-auto w-full m-1 box-border">
      <MonacoEditor
        v-model:value="strText"
        :options="options"
        :height="editorH - 20"
        :weight="editorW"
      />
    </div>
  </div>
</template>

<style>
</style>
