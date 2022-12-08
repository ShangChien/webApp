<script setup lang="ts">
import { defineComponent,ref,h, } from 'vue'
import type {CSSProperties} from 'vue'
import { repeat } from 'seemly'
import gridPage from "@/components/rdkitComponent/gridPage.vue";
import type { TreeOption } from 'naive-ui'
import { NTree,NButtonGroup,NButton,NIcon,NSwitch,NSpace } from 'naive-ui'
import { useEnumStore } from '@/stores/enumStore'
const enumStore = useEnumStore();
const props = defineProps();
const emit = defineEmits([]);
function createData (level = 4, baseKey = ''): TreeOption[] | undefined {
  if (!level) return undefined
  return repeat(6 - level, undefined).map((_, index) => {
    const key = '' + baseKey + level + index
    return {
      label: createLabel(level),
      key,
      children: createData(level - 1, key)
    }
  })
}

function createLabel (level: number): string {
  if (level === 4) return '道生一'
  if (level === 3) return '一生二'
  if (level === 2) return '二生三'
  if (level === 1) return '三生万物'
  return ''
}
const data=createData()
const defaultExpandedKeys = ref(['40', '41'])
const renderSwitcherIconWithExpaned=({ expanded }: { expanded: boolean }) =>
        h(NIcon, null, {
          default: () => h(expanded ? h('div',{class:'i-fxemoji-openfolder ml--0.9 mt--1.5 rotate--90 text-xl'}) : h('div',{class:'i-fxemoji-folder text-xl mt--1'}))
        })
const checkStrategy = ref<'label'|'type'>('label')
const railStyle=({focused,checked}:{focused:boolean;checked:boolean})=>{
  const style: CSSProperties = {}
  if (checked) {
    style.background = '#7dd3fc'
    if (focused) {
      style.boxShadow = '0 0 0 2px #7dd3fc40'
    }
  } else {
    style.background = '#6ee7b7'
    if (focused) {
      style.boxShadow = '0 0 0 2px #6ee7b740'
    }
  }
  return style
}
const fun=(value:string)=>console.log(value)
</script>

<template>
<div class="b-2 rd-2 b-indigo-100">
<n-space justify="space-between">
  <n-switch size="large"
            :rail-style="railStyle" 
            :round="false"
            checked-value="label"
            unchecked-value="type"
            v-model:value=checkStrategy
            @click="fun(checkStrategy)" >
    <template #checked>group label</template>
    <template #checked-icon><div class="i-fluent-emoji-flat-label text-2xl" /></template>
    <template #unchecked>group type</template>
    <template #unchecked-icon><div class="i-fluent-emoji-bubbles text-xl" /></template>
  </n-switch>
  <n-space>
  <n-button-group size='small'>
    <n-button value="label">
      fun1
    </n-button>
    <n-button value="type">
      fun2
    </n-button>
  </n-button-group>
  </n-space>
</n-space>
<n-tree 
	block-line
  :data="data"
  :default-expanded-keys="defaultExpandedKeys"
  :render-switcher-icon="renderSwitcherIconWithExpaned"
  selectable expand-on-click
/>
</div>
</template>
<style>
</style>