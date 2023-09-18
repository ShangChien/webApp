<script setup lang='ts'>
import { ref, watch } from 'vue'
import { computedEager, useManualRefHistory } from '@vueuse/core'
import monaco from '@/components/monaco/monaco.vue'
import editorHistory from '@/components/monaco/editorHistory.vue'

// change strText.value to update children text
// childrenDom inner change can be watched by text4View
// records.history.value[0].snapshot.text is latest textstr
// records.history.value.slice(-1)[0].snapshot.text(= props.strText; strText.value) is init textstr
const props = defineProps<{ strText: string; name: string }>()
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
  records.commit()
  detail.value = ''
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
    <div class="flex-none w-full flex flex-nowrap justify-between items-center pt-1 px-1 box-border gap-2">
      <div class="flex-auto flex flex-nowrap justify-start items-center box-border gap-2">
        <div
          v-if="canSave"
          class="bg-sky-2 rd-1 p-1 m-0 text-(x blue-5) leading-1em cursor-pointer box-border
              hover:(bg-sky-3)
              active:(outline outline-2px outline-blue-3)"
          @click="save()"
        >
          save
        </div>
        <div
          v-if="isDiff"
          class="bg-sky-2 rd-1 p-1 m-0 text-(x blue-5) leading-1em cursor-pointer box-border
              hover:(bg-sky-3)
              active:(outline outline-2px outline-blue-3)"
          @click="commit()"
        >
          commit
        </div>
        <div
          v-if="records.canUndo.value"
          class="bg-sky-2 rd-1 p-1 m-0 text-(x blue-5) leading-1em cursor-pointer box-border
              hover:(bg-sky-3)
              active:(outline outline-2px outline-blue-3)"
          @click="undo()"
        >
          undo
        </div>
        <div
          v-if="records.canRedo.value"
          class="bg-sky-2 rd-1 p-1 m-0 text-(x blue-5) leading-1em cursor-pointer box-border
              hover:(bg-sky-3)
              active:(outline outline-2px outline-blue-3)"
          @click="redo()"
        >
          redo
        </div>
        <div
          v-if="canReset"
          class="bg-sky-2 rd-1 p-1 m-0 text-(x blue-5) leading-1em cursor-pointer box-border
              hover:(bg-sky-3)
              active:(outline outline-2px outline-blue-3)"
          @click="reset()"
        >
          reset
        </div>
      </div>
      <div class="flex-none leading-1em">
        {{ props.name ? props.name : "empty (o_o)" }}
      </div>
      <editorHistory
        :records="records.history.value"
        @backTo="(index) => console.log('back to', index)"
        @diffTo="(index) => console.log('diff to', index)"
      />
    </div>
    <monaco ref="monacoDom" :key="RefreshKey" v-model:currentText="text4View.text" :strText="strText" />
  </div>
</template>

<style>
</style>
