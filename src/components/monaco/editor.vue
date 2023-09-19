<script setup lang='ts'>
import { ref, watch } from 'vue'
import { computedEager, useManualRefHistory } from '@vueuse/core'
import { NButton } from 'naive-ui'
import monaco from '@/components/monaco/monaco.vue'
import editorHistory from '@/components/monaco/editorHistory.vue'

// change strText.value to update children text
// childrenDom inner change can be watched by text4View
// records.history.value[0].snapshot.text is latest textstr
// records.history.value.slice(-1)[0].snapshot.text(= props.strText; strText.value) is init textstr
const props = defineProps<{ strText: string; name: string }>()
const emits = defineEmits<{ sync: [strText: string] }>()
const monacoDom = ref()
const RefreshKey = ref<number>(0)
const strText = ref<string>(props.strText)
const detail = ref<string>('')
const text4View = ref<{ text: string; detial: string }>({
  text: strText.value,
  detial: detail.value,
})
const records = useManualRefHistory(text4View, { clone: true })
const isDiff = computedEager(() => {
  return records.history.value[0].snapshot.text !== text4View.value.text
})
const canSave = computedEager(() => {
  return (props.strText !== strText.value) || (props.strText !== text4View.value.text)
})
const canReset = computedEager(() => {
  return records.history.value.slice(-1)[0]?.snapshot.text !== text4View.value.text
})

function save() {
  text4View.value.detial = `save: ${detail.value}`
  detail.value = ''
  records.commit()
  emits('sync', text4View.value.text)
}

function commit() {
  text4View.value.detial = `commit: ${detail.value}`
  records.commit()
  detail.value = ''
}

function undo() {
  records.undo()
  strText.value = text4View.value.text
}

function redo() {
  records.redo()
  strText.value = text4View.value.text
}
function reset() {
  text4View.value.text = props.strText
  strText.value = props.strText
  records.commit()
  records.clear()
  RefreshKey.value++
}

watch(() => props.strText, (val) => {
  // do other data operation before records.clear(), if necessary
  strText.value = val
  text4View.value = { text: val, detial: '' }
  detail.value = ''
  records.commit()
  records.clear()
})
</script>

<template>
  <div class="flex flex-col flex-nowrap justify-between items-center h-full w-full box-border p-1">
    <div class="flex-none w-full flex flex-nowrap justify-between items-center box-border gap-0">
      <div class="flex-auto flex flex-nowrap justify-start items-center box-border gap-2 p-0">
        <NButton size="small" :type="!canSave ? 'success' : 'warning'" :disabled="!canSave" @click="save()">
          {{ !canSave ? "data safe" : "need sync" }}
        </NButton>
        <NButton size="small" type="info" :disabled="!isDiff" @click="commit()">
          local commit
        </NButton>
        <NButton size="small" secondary type="info" :disabled="!records.canUndo.value" @click="undo()">
          undo
        </NButton>
        <NButton size="small" secondary type="info" :disabled="!records.canRedo.value" @click="redo()">
          redo
        </NButton>
        <NButton size="small" type="error" :disabled="!canReset" @click="reset()">
          reset
        </NButton>
      </div>
      <div class="flex-none leading-1em mx-1">
        {{ props.name ? props.name : "" }}
      </div>
      <editorHistory
        :records="records.history.value"
        @backTo="(index) => console.log('back to', index)"
        @diffTo="(index) => console.log('diff to', index)"
      />
    </div>
    <div class="w-full flex-auto editorH br-1 box-border">
      <monaco ref="monacoDom" :key="RefreshKey" v-model:currentText="text4View.text" :strText="strText" />
    </div>
  </div>
</template>

<style>
.editorH {
  height: calc(100% - 36px)
}
</style>
