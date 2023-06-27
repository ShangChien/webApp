<script setup lang="ts">
import { computed, h, inject, onMounted, ref } from 'vue'
import type { Component, Ref } from 'vue'
import {
  NButton,
  NCard,
  NCheckbox,
  NDropdown,
  NIcon,
  NPopover,
  NTag,
} from 'naive-ui'
import { useClipboard } from '@vueuse/core'
import { Dots } from '@vicons/tabler'
import { CopyFile, Delete, Edit } from '@vicons/carbon'
import svgRdkit from '@/components/rdkitComponent/svgRdkit.vue'
import { useCopy } from '@/components/rdkitComponent/composable/useCopy'
import type { pgDataItem } from '@/components/types'
import { useEnumStore } from '@/stores/enumStore'

const props = defineProps<pgDataItem>()
const _emit = defineEmits(['itemChecked'])
const currentEdit: Ref<{ id: number;state: number }> = inject('currentEdit')
const enumStore = useEnumStore()
const editState = computed(() => {
  if (currentEdit.value.state === 2 && currentEdit.value.id === props.id)
    return 'editing'

  else if (currentEdit.value.state === 1 && currentEdit.value.id === props.id)
    return 'preview'

  else
    return 'normal'
})

// 可视加载组件
const checked = ref(false)
const cardView = ref()
const target = ref(null)

const { copy } = useClipboard()
const _copytext = ref<string>(props.smiles)

function renderIcon(icon: Component) {
  return () => {
    return h(NIcon, null, {
      default: () => h(icon),
    })
  }
}
const options = [
  { label: 'edit', key: 'edit', icon: renderIcon(Edit) },
  { label: 'copy', key: 'copy', icon: renderIcon(CopyFile) },
  { label: 'delete', key: 'delete', icon: renderIcon(Delete) },
]
function handleDropOption(key: string) {
  console.log(key)
  switch (key) {
    case 'edit':
      currentEdit.value.id = props.id
      currentEdit.value.state = 1
      // console.log('editing',editState.value)
      break
    case 'copy':
      useCopy('')
      // copy(copytext.value);
      console.log('copied')
      break
    case 'delete':
      enumStore.rmMolById(props.id)
      // console.log('delete ',props.id)
      break
    default:
      console.log('no option matched ')
  }
}
function quickView() {
  if (currentEdit.value.id !== props.id) {
    currentEdit.value.id = props.id
    currentEdit.value.state = 1
    console.log('preview', editState.value)
  }
  else {
    currentEdit.value.id = 0
    currentEdit.value.state = 0
    console.log('exit preview')
  }
}
onMounted(() => {
  // console.log(props);
})
</script>

<template>
  <div ref="target" class="w-100% h-100% ">
    <NCard ref="cardView" hoverable>
      <template #cover>
        <div class="ml-1 z-1 w-100%  flex justify-between">
          <NCheckbox v-model:checked="checked" />
          <div>
            <NTag v-if="editState === 'editing'" :bordered="false" type="warning" size="tiny">
              {{ `editing:${props.id}` }}
            </NTag>
            <NTag v-else-if="editState === 'preview'" :bordered="false" type="primary" size="tiny">
              {{ `preview:${props.id}` }}
            </NTag>
            <NTag v-else :bordered="false" type="info" size="tiny">
              {{ `id: ${props.id}` }}
            </NTag>
          </div>
          <NDropdown
            trigger="click"
            :options="options"
            placement="bottom-end"
            @select="handleDropOption"
          >
            <NButton tertiary circle size="tiny" class="scale-90 mr-1">
              <template #icon>
                <NIcon><Dots /></NIcon>
              </template>
            </NButton>
          </NDropdown>
        </div>
        <NPopover
          trigger="hover"
          placement="right"
          display-directive="if"
          to=".n-scrollbar"
        >
          <template #trigger>
            <div class="w-100% h-100%" @dblclick="quickView">
              <svg-rdkit v-bind="props" />
            </div>
          </template>
          <span
            class="cursor-pointer rd-1 hover:bg-blue-50"
            @click="copy(props.name_calc)"
          >计算:{{ props.name_calc }}</span><br>
          <span
            class="cursor-pointer rd-1 hover:bg-blue-50"
            @click="copy(props.name_mat)"
          >材料:{{ props.name_mat }}</span>
        </NPopover>
      </template>
    </NCard>
  </div>
</template>

<style></style>
