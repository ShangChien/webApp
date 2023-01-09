<script setup lang="ts">
import { inject,ref,h,computed,watch } from 'vue'
import type {CSSProperties,VNode,Ref,ComputedRef} from 'vue'
import gridPage from "@/components/rdkitComponent/gridPage.vue";
import type { TreeOption } from 'naive-ui'
import { NTree,NSwitch,NSpace } from 'naive-ui'
import { Splitpanes, Pane } from 'splitpanes'
import { useMouseInElement } from '@vueuse/core'
import { useEnumStore } from '@/stores/enumStore'
import type { molData } from "@/components/types"
const enumStore = useEnumStore();
const siderBox=ref(null)
const { isOutside } = useMouseInElement(siderBox)
const checkStrategy = ref<'label'|'type'>('type')
const expandedKeys = ref([])
const checkedKeys = ref([])
function unique(arr:any[],val:string) {
    const res = new Map();
    return arr.filter(item => item!=undefined ? (!res.has(item[val]) && res.set(item[val],1)) : false)
  }
const gridData = computed(()=>{
    let out:molData[] = checkedKeys.value.map((v) => v.includes('-') ?
      enumStore.getById(+v.split('-').at(-1))[0] : undefined
    )
    return unique(out,'id')
})
const currentEdit:Ref<{id:number;state:number}>=inject('currentEdit')
//refresh checkedKey
watch(
  ()=>currentEdit.value.state,
  ()=>{
    if (currentEdit.value.state==3) {
      console.log('watch currentid-1',currentEdit.value.id,currentEdit.value.state,checkedKeys.value)
      checkedKeys.value = checkedKeys.value.filter((v) => v.includes('-'))
      console.log('watch currentid-2',currentEdit.value.id,currentEdit.value.state,checkedKeys.value)
      currentEdit.value.state=0
    }
  }
)
const labels:ComputedRef<string[]>=computed(()=>enumStore.getAllLabels)
//const labels=enumStore.getAllLabels
const types=['mole','ligand','core']
const siderData =computed<TreeOption[]|undefined>(()=>{
  if (checkStrategy.value=='label') {
    return labels.value.map((label,index)=>{
      const key = ''+index
      const children= enumStore.getByLabel(label).map((mol,indexInner)=>{
        return{
          label : ((mol?.name)??(''+mol.id))+'-'+mol.id,
          key : key+'-'+indexInner+'-'+mol.id,
        }
      })
      return { label,key,children }
    })
  } else {
    return types.map((type:'core'|'ligand'|'mole',index)=>{
      const key = ''+index
      const children= enumStore.getByType(type).map((mol,indexInner)=>{
        return{
          label : ((mol?.name)??(''+mol.id))+'-'+mol.id,
          key  : key+'-'+indexInner+'-'+mol.id,
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
  expandedKeys.value = []
  checkedKeys.value = []
  currentEdit.value.id=0
  currentEdit.value.state=0
}
const onCheck = (key:string[])=>{
  checkedKeys.value=key
  //console.log(checkedKeys.value)
}
const onSelect = (key:string[])=>{
  currentEdit.value.state=0
  let v = key[0]
  if ((key[0]?.includes('-'))) {
    currentEdit.value.id= +v.split('-').at(-1)
    currentEdit.value.state=1
  }else{
    currentEdit.value.id=0
  }
  console.log('currentEdit:',currentEdit.value)
}
defineExpose({gridData})
</script>

<template>
<splitpanes class="pt-1 splitpanes h-full " style="background-color: #ffffff">
  <pane size="15" min-size="11">
    <div ref="siderBox" class="b-2 rd-2 b-indigo-100 min-w-180px h-89vh">
      <n-space justify="space-between" class="p-1" >
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
          <div class="text-center rd-1 w-5 hover:(bg-gray-200 cursor-pointer)" @click="">
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
        :checked-keys="checkedKeys"
        :render-prefix="renderPrefix"
        :render-switcher-icon="renderSwitcherIcon"
        :node-props="nodeProps"
        @update:checked-keys="onCheck"
        @update:selected-keys="onSelect"
        cascade checkable virtual-scroll />
    </div>
  </pane>
  <pane size="60" min-size="16" >
    <grid-page v-if="gridData" :molList="gridData" :cols='6' :rows="6" class="h-89vh" />
  </pane>
  <pane size="25" min-size="10" >
    <div class="w-full h-full text-center">
      <div class="v-mid leading-700px" >
        <div class="i-vscode-icons-file-type-flareact text-5xl"></div>
      </div>
    </div>
  </pane>
</splitpanes>
</template>
<style >
@import "splitpanes/dist/splitpanes.css";
.splitpanes {background-color: #f8f8f8;}

.splitpanes__splitter {
  background-color:#ecfeff;
  position: relative;
  width:3px;
  border-radius: 0.2rem;
  margin: 0.2rem;
}
.splitpanes__splitter:before {
  content: '';
  position: absolute;
  left: 0;
  top: 0;
  border-radius: 0.4rem;
  transition: opacity 0.4s;
  background-color: #6bd8fca3;
  opacity: 0;
  z-index: 1;
}
.splitpanes__splitter:hover:before {opacity: 1;}
.splitpanes--vertical > .splitpanes__splitter:before {left: -2px;right: -2px;height: 100%;}
.splitpanes--horizontal > .splitpanes__splitter:before {top: -2px;bottom: -2px;width: 100%;}
</style>