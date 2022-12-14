<script setup lang="ts">
import { ref,h,computed,nextTick,toRaw } from 'vue'
import type {CSSProperties,VNode} from 'vue'
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
const expandedKeys = ref([])
const checkedKeys = ref([])
const gridData=ref([])
const labels:string[] = enumStore.getAllLabels
const types=['mole','ligand','core']
const siderData =computed<TreeOption[]|undefined>(()=>{
  if (checkStrategy.value=='label'){
    return labels.map((label,index)=>{
      const key = ''+index
      const children= enumStore.getByLabel(label).map((mol,indexInner)=>{
        return{
          label : (mol?.name)??(''+mol.id),
          key : key+'-'+indexInner,
        }
      })
      return { label,key,children }
    })
  } else {
    return types.map((type:'core'|'ligand'|'mole',index)=>{
      const key = ''+index
      const children= enumStore.getByType(type).map((mol,indexInner)=>{
        return{
          label : (mol?.name)??(''+mol.id),
          key  : key+'-'+indexInner,
        }
      })
      return { label:type,key,children }
    })
  }
})
const renderPrefix=({ option }: { option: TreeOption }):VNode=>{
  let render:VNode
  if ((option.key as string).includes('-')) {
    render = h('div',{ class:'i-flat-color-icons-mind-map' })
  }else if (expandedKeys.value.includes(option.key)) {
    render = h('div',{class:'i-fxemoji-openfolder text-xl'}) 
  }else{
    render = h('div',{class:'i-fxemoji-folder text-xl '})
  }
  return render
}
const nodeProps=({ option }: { option: TreeOption }) => {
        return {
          onDblclick () {
            if (!(option.key as string).includes('-')) {
            expandedKeys.value.includes(option.key)?
            expandedKeys.value.splice(expandedKeys.value.indexOf(option.key),1):
            expandedKeys.value.push(option.key)
            console.log( expandedKeys.value)}
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
const switchClass=()=>{
  expandedKeys.value=[]
}
const onCheck=(key:string[])=>{
  checkedKeys.value=key
  console.log('checkkey:',key)
}

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
              @click="switchClass" >
      <template #checked>types</template>
      <template #checked-icon><div class="i-fluent-emoji-bubbles text-xl" /></template>
      <template #unchecked>labels</template>
      <template #unchecked-icon><div class="i-fluent-emoji-flat-label text-2xl" /></template>
    </n-switch>
    <div v-show="!isOutside" class="grid grid-cols-4 ">
      <div class="text-center rd-1 w-5 hover:(bg-gray-200 cursor-pointer)">
        <div class="i-codicon-new-file text-4"></div>
      </div>
      <div class="text-center rd-1 w-5 hover:(bg-gray-200 cursor-pointer)">
         <div class="i-codicon-new-folder text-4"></div>
      </div>
      <div class="text-center rd-1 w-5 hover:(bg-gray-200 cursor-pointer)">
        <div class="i-codicon-refresh text-4"></div>
      </div>
      <div class="text-center rd-1 w-5 hover:(bg-gray-200 cursor-pointer)" @click="expandedKeys=[]">
        <div class="i-codicon-collapse-all text-4"></div>
      </div>
    </div>
  </n-space>
  <n-tree class="h-86vh"
  	:block-line='true' :block-node='true'
    :data="siderData"
    :expanded-keys="expandedKeys"
    :render-prefix="renderPrefix"
    :render-switcher-icon="renderSwitcherIcon"
    :node-props="nodeProps"
    @update:checked-keys="onCheck"
    cascade checkable virtual-scroll
     />
  </div>
  <grid-page v-if="gridData" :molList="gridData" :cols='6' class="h-89.8vh" />
</div>
</template>
<style>
.gridlayout {
  display: grid;
  grid-column-gap: 0.5em;
  grid-template-columns: 2fr 6fr 4fr;
}
</style>