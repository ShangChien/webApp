<script setup lang="ts">
import { defineAsyncComponent, h, onMounted, reactive, ref } from 'vue'
import {
  NButton,
  NCard,
  NCheckbox,
  NDropdown,
  NEmpty,
  NIcon,
  NModal,
  NPopover,
  NSpace,
  NSpin,
} from 'naive-ui'
import { useClipboard, useMouseInElement } from '@vueuse/core'
import type { molData } from '@/components/types'
import svgRdkit from '@/components/rdkitComponent/svgRdkit.vue'

const props = defineProps<molData>()
const emit = defineEmits(['update:smiles',
  'update:atoms',
  'update:bonds',
  'update:selected',
  'update:label',
  'itemDeleted'])
// 异步加载mol编辑组件
const editableRdkit = defineAsyncComponent({
  loader: () => import('@/components/rdkitComponent/editableRdkit.vue'),
  loadingComponent: NSpin,
  delay: 2000, // +600*Math.random(),
  errorComponent: NEmpty,
  timeout: 3000,
  suspensible: false,
})
// const { smiles,atoms,bonds,selected,label } = useVModels(props, emit)
// 可视加载组件
const checked = ref(false)
const cardView = ref()
const { isOutside } = useMouseInElement(cardView)
const target = ref(null)

// handle Modal
const showModal = ref(false)
function onNegativeClick() {
  showModal.value = false
}

function onPositiveClick() {
  showModal.value = false
}
const initSmile = reactive({
  smiles: props.smiles,
  atoms: ref(props.atoms),
  bonds: ref(props.bonds),
})
function acceptMol(mol: molData) {
  initSmile.smiles = mol.smiles
  initSmile.atoms = mol.atoms
  initSmile.bonds = mol.bonds
  console.log(initSmile, 'cardRdkit.vue acceptMol')
}

const { copy } = useClipboard()
const copytext = ref<string>(props.smiles)

function renderIcon(iconClassStr: string) {
  return () => h(NIcon, null, { default: () => h('div', { class: iconClassStr }) })
}
const options = [
  { label: 'edit', key: 'edit', icon: renderIcon('i-carbon-edit') },
  { label: 'copy', key: 'copy', icon: renderIcon('i-carbon-copy-file') },
  { label: 'delete', key: 'delete', icon: renderIcon('i-carbon-trash-can') },
]
function handleDropOption(key: string) {
  if (key === 'edit') {
    showModal.value = true
  } else if (key === 'copy') {
    copy(copytext.value)
  } else if (key === 'delete') {
    emit('itemDeleted')
  } else if (key === 'edit') {
    //
  }
}
onMounted(() => {
  console.log(props)
})
</script>

<template>
  <div ref="target" style="width:100%;height:100%">
    <!-- v-if="targetIsVisible" -->
    <NCard
      ref="cardView"
      hoverable
    >
      <template #cover>
        <NSpace
          v-if="!isOutside || checked"
          justify="space-between"
          style="
          padding-left: 3%;
          padding-top: 2%;
          padding-right: 1%;
          position: absolute;
          z-index: 1;
          width: 95%;
        "
        >
          <NCheckbox v-model:checked="checked" style="opacity: 0.8;" />
          <div>
            <NDropdown
              trigger="click"
              :options="options"
              placement="bottom-end"
              @select="handleDropOption"
            >
              <NButton tertiary circle size="tiny">
                <template #icon>
                  <NIcon><div class="i-tabler-dots" /></NIcon>
                </template>
              </NButton>
            </NDropdown>
          </div>
        </NSpace>
        <NModal v-model:show="showModal" :mask-closable="false">
          <NCard
            style="width: 900px; height: 800px"
            title="位点重新选取"
            :bordered="false"
            size="huge"
            role="dialog"
            aria-modal="true"
          >
            <editable-rdkit
              v-bind="props"
              qsmiles="*~*"
              style="width: 70%"
              @update-mol="acceptMol"
            />
            <template #footer>
              <NSpace>
                <NButton @click="onNegativeClick">
                  取消
                </NButton>
                <NButton @click="onPositiveClick">
                  确定
                </NButton>
              </NSpace>
            </template>
          </NCard>
        </NModal>
        <NPopover
          trigger="hover"
          placement="right"
          width="trigger"
          display-directive="if"
          to=".n-scrollbar"
        >
          <template #trigger>
            <div style="width: 100%; height: 100%">
              <svg-rdkit v-bind="props" />
            </div>
          </template>
          <span style="word-break: break-word">
            {{ props.smiles }}<br>
            atoms:{{ props.atoms }}<br>
            bonds:{{ props.bonds }}
          </span>
        </NPopover>
      </template>
    </NCard>
  </div>
</template>

<style></style>
