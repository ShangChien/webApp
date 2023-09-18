<!-- eslint-disable new-cap -->
<script setup lang='ts'>
import { onMounted, ref, watch } from 'vue'
import * as monaco from 'monaco-editor'
import editorWorker from 'monaco-editor/esm/vs/editor/editor.worker?worker'
import jsonWorker from 'monaco-editor/esm/vs/language/json/json.worker?worker'
import cssWorker from 'monaco-editor/esm/vs/language/css/css.worker?worker'
import htmlWorker from 'monaco-editor/esm/vs/language/html/html.worker?worker'
import tsWorker from 'monaco-editor/esm/vs/language/typescript/ts.worker?worker'
import lightTheme from '@/assets/light-vitesse.json'
import softLightTheme from '@/assets/soft-light-vitesse.json'

const props = defineProps<{ strText: string }>()
const currentText = defineModel<string>('currentText')
const editorDom = ref(null)

// eslint-disable-next-line no-restricted-globals
self.MonacoEnvironment = {
  getWorker(_, label) {
    if (label === 'json') {
      return new jsonWorker()
    }
    if (label === 'css' || label === 'scss' || label === 'less') {
      return new cssWorker()
    }
    if (label === 'html' || label === 'handlebars' || label === 'razor') {
      return new htmlWorker()
    }
    if (label === 'typescript' || label === 'javascript') {
      return new tsWorker()
    }
    return new editorWorker()
  },
}

onMounted(() => {
  // @ts-expect-error: set theme.json
  monaco.editor.defineTheme('vitesse-light', lightTheme)
  // @ts-expect-error: set theme.json
  monaco.editor.defineTheme('soft-vitesse-light', softLightTheme)
  const editor = monaco.editor.create(editorDom.value, {
    colorDecorators: true,
    lineHeight: 1.2,
    tabSize: 2,
    theme: 'vitesse-light',
    language: 'python',
    automaticLayout: true,
  })

  // updata editor contents by
  watch(() => props.strText, (val) => {
    editor.setValue(val)
    // editor.updateOptions({ theme: 'vs' })
  }, { immediate: true })

  editor.onDidChangeModelContent((_event) => {
    currentText.value = editor.getValue()
  })
})
</script>

<template>
  <div ref="editorDom" class="h-full w-full p-1 box-border" />
</template>

<style>
</style>
