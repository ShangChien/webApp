<script setup lang="ts">
import { ref,h,computed } from 'vue'
import type {CSSProperties} from 'vue'
import gridPage from "@/components/rdkitComponent/gridPage.vue";
import type { TreeOption } from 'naive-ui'
import { NTree,NSwitch,NSpace } from 'naive-ui'
import { useMouseInElement } from '@vueuse/core'
import { useEnumStore } from '@/stores/enumStore'
import type { molData } from "@/components/types"
const enumStore = useEnumStore();
const props = defineProps();
const emit = defineEmits([]);
const siderBox=ref(null)
const { isOutside } = useMouseInElement(siderBox)
const checkStrategy = ref<'label'|'type'>('label')
const defExpandedKeys = ref([])
const labels:string[] = enumStore.getAllLabels
const types=['mole','ligand','core']
const useData =computed<TreeOption[]|undefined>(()=>{
  if (checkStrategy.value=='label'){
    return labels.map((label,index)=>{
      const key = index
      const prefix =()=> h(defExpandedKeys.value.includes(key) ? 
        h('div',{class:'i-fxemoji-openfolder text-xl'}) 
        :h('div',{class:'i-fxemoji-folder text-xl '}))
      const children= enumStore.getByLabel(label).map((mol,indexInner)=>{
        return{
          label : (mol?.name)??(''+mol.id),
          key : ''+index+indexInner,
          prefix : () => h('div',{ class:'i-flat-color-icons-mind-map' })
        }
      })
      return { label,key,prefix,children }
    })
  } else {
    return types.map((type:'core'|'ligand'|'mole',index)=>{
      const key = index
      const prefix =()=> h(defExpandedKeys.value.includes(key) ? 
        h('div',{class:'i-fxemoji-openfolder text-xl'}) 
        :h('div',{class:'i-fxemoji-folder text-xl '}))
      const children= enumStore.getByType(type).map((mol,indexInner)=>{
        return{
          label : (mol?.name)??(''+mol.id),
          key  :''+index+indexInner,
          prefix : () => h('div',{ class:'i-flat-color-icons-mind-map' })
        }
      })
      return { label:type,key,prefix,children }
    })
  }
})
const nodeProps=({ option }: { option: TreeOption }) => {
        return {
          onClick () {
            defExpandedKeys.value.includes(option.key)?
            defExpandedKeys.value.splice(defExpandedKeys.value.indexOf(option.key),1):
            defExpandedKeys.value.push(option.key)
          },
          onContextmenu (e: MouseEvent): void {
            e.preventDefault()
          }
        }
      }
const renderSwitcherIcon=()=>h('div', { class:'i-ion-chevron-forward'})
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
const sort=(value:string)=>{
  console.log()
}

const ligands=ref<molData[]|any>()
const cores=ref<molData[]|any>()
ligands.value=enumStore.getByType('ligand')
cores.value=enumStore.getByType('core')
</script>

<template>
<div class="pt-1 gridlayout">
  <div ref="siderBox" class="b-2 rd-2 p-1 b-indigo-100 min-w-170px" >
  <n-space justify="space-between" >
    <n-switch size="medium"
              :rail-style="railStyle" 
              :round="false"
              checked-value="type"
              unchecked-value="label"
              v-model:value=checkStrategy
              @click="sort(checkStrategy)" >
      <template #checked>types</template>
      <template #checked-icon><div class="i-fluent-emoji-bubbles text-xl" /></template>
      <template #unchecked>labels</template>
      <template #unchecked-icon><div class="i-fluent-emoji-flat-label text-2xl" /></template>
    </n-switch>
    <div v-show="!isOutside" class="grid grid-cols-4 ">
      <div class="text-center rd-1 w-5 hover:bg-gray-200">
        <div class="i-codicon-new-file text-4"></div>
      </div>
      <div class="text-center rd-1 w-5 hover:bg-gray-200">
         <div class="i-codicon-new-folder text-4"></div>
      </div>
      <div class="text-center rd-1 w-5 hover:bg-gray-200">
        <div class="i-codicon-refresh text-4"></div>
      </div>
      <div class="text-center rd-1 w-5 hover:bg-gray-200">
        <div class="i-codicon-collapse-all text-4"></div>
      </div>
    </div>
  </n-space>
  <n-tree 
  	:block-line='true' :block-node='true'
    :data="useData"
    :default-expanded-keys="defExpandedKeys"
    :render-switcher-icon="renderSwitcherIcon"
    :node-props="nodeProps"
    selectable />
  </div>
  <grid-page :molList="ligands" :cols='6' class="h-89.8vh" />
</div>
</template>
<style>
.gridlayout {
  display: grid;
  grid-column-gap: 0.5em;
  grid-template-columns: 2fr 6fr 4fr;
}
</style>