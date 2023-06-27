<script setup lang="ts">
import { NSelect, NTag } from 'naive-ui'
import { computed, h, onMounted } from 'vue'
import { useVModel } from '@vueuse/core'
import { useEnumStore } from '@/stores/enumStore'

const props = defineProps<{ tags: string[] }>()
const emit = defineEmits(['update:tags'])
const tagsValue = useVModel(props, 'tags', emit)

const enumStore = useEnumStore()
const tagsOption = computed<{ label: string; value: string }[]>(() => {
  return enumStore.getAllLabels.map((v: string) => ({ label: v, value: v }))
})
function renderTag({ option, handleClose }) {
  return h(NTag, {
    type: 'success',
    closable: true,
    round: true,
    size: 'small',
    onMousedown: (e: FocusEvent) => {
      e.preventDefault()
    },
    onClose: (e: MouseEvent) => {
      e.stopPropagation()
      handleClose()
    },
  },
  { default: () => option.label })
}
onMounted(() => {

})
</script>

<template>
  <NSelect
    v-model:value="tagsValue"
    type="info"
    size="small"
    :options="tagsOption"
    :render-tag="renderTag"
    multiple clearable filterable tag
  />
</template>
