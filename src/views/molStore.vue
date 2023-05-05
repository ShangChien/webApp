<script setup lang="ts">
import { computed,ref,reactive } from 'vue'
import { NButton,NCheckbox,NSwitch } from 'naive-ui'
import { createReusableTemplate } from '@vueuse/core'
import { useMolStore } from '@/stores/molStore'
import type { molData } from "@/components/types"
import gridPage from "@/components/rdkitComponent/gridPage.vue";
interface dataPick {
  items: string[],
  selected: string[],
  expend:boolean
  enabled:boolean
}
const angle=ref(0)
const [DefineTag, ReuseTag] = createReusableTemplate<{label:string }>()
const store = useMolStore()
const labels = computed(()=>store.getAllLabels)
const dataLabels:dataPick=reactive({items:labels,selected:[],expend:false,enabled:false})
const dataTypes:dataPick=reactive({items:['molecule','ligand','core'],selected:[],expend:false,enabled:false})
const mols4view = computed(()=>{
  return store.getByTypesAndLabels(dataTypes.selected,dataLabels.selected)
})
async function selection(data:dataPick) {
  let len = await Promise.resolve(data.items.length)
  if (len !== 0){
    data.selected=(data.selected.length!==0) ? [] : [...data.items]
  }else{
    console.log('data items length', data.selected.length)
  }
}
async function reverse(data:dataPick) {
  //翻转icon
  angle.value = (angle.value == 0) ? 180 : 0
  let arr= await Promise.resolve(data.items)
  data.selected = arr.filter(item => !data.selected.includes(item))
}

</script>

<template>
<div>
  <define-tag v-slot="{label}">
    <n-button class="m-0.5" round strong secondary size='small'
      :type="dataLabels.selected.includes(label) ? 'primary':'tertiary'"
      @click="dataLabels.selected.includes(label) ? 
      dataLabels.selected.splice(dataLabels.selected.indexOf(label), 1) : dataLabels.selected.push(label)">
      {{label}}
    </n-button>
  </define-tag>
  <div class="filter">
    <div class="flex flex-nowrap justify-between items-start">
      <div class="flex-none flex items-center mr-2 rd-2 pr-1 pl-1 mb-0.5 mt-0.5 bg-blue-100 cursor-pointer box-border b-2 transition-210"
      :class="[dataLabels.enabled 
      ? 'bg-indigo-400 b-indigo-300 hover:(bg-indigo-300 b-indigo-400)'
      : 'bg-indigo-300 b-indigo-400 hover:(bg-indigo-400 b-indigo-300)']"
      @click="dataLabels.enabled=!dataLabels.enabled">
        <div class="text-2xl mt--0.5 mb-0.5"
        :class="[dataLabels.enabled ? 'i-fluent-emoji-flat-label' : 'i-fluent-emoji-flat-label?mask text-zinc-100']"></div>
        <span class="text-center text-zinc-100">标签</span>
      </div>
      <Transition name="fade">
        <div v-if="dataLabels.enabled" class="flex flex-nowrap justify-begin overflow-hidden flex-auto items-start">
          <div class="rd-2 mt-0.5 bg-blue-100 flex-none flex flex-nowrap justify-around items-center hover:bg-blue-200 ">
            <n-checkbox class=" ml-1.5 mr-1" size="small"
            :checked="dataLabels.selected.length!==0 && dataLabels.selected.length===dataLabels.items.length"
            :indeterminate="dataLabels.selected.length!==0 && dataLabels.selected.length!==dataLabels.items.length"
            @click="selection(dataLabels)" >
            </n-checkbox>
            <div class="i-carbon-contrast rd-2 m-1 ml-2  bg-indigo-400 cursor-pointer text-xl transition-210 hover:(bg-indigo-600)"
            :style="{ 'transform': `rotateY(${angle}deg)` }"
            @click="reverse(dataLabels)"></div>
          </div>  
          <div class="flex justify-begin overflow-hidden flex-auto"
            :class="[dataLabels.expend ? 'flex-wrap':'flex-nowrap']">
            <reuse-tag v-for="label in dataLabels.items" :label="label"></reuse-tag>
          </div>
          <div class="rd-2 ml-2 p-0.5 mt-0.5 bg-blue-100 flex-none hover:bg-blue-200">
            <div @click="dataLabels.expend = !dataLabels.expend"
            :class="[dataLabels.expend ? 'i-ion-chevron-collapse' : 'i-ion-chevron-expand']"
            class=" text-l p-1 c-zinc-400 hover:(c-slate-600 cursor-pointer)" />
          </div>
        </div>
      </Transition>
    </div>
    <div class="flex flex-nowrap justify-between items-start">
      <div class="flex-none flex items-center mr-2 rd-2 pr-1 pl-1 mb-0.5 mt-0.5 bg-blue-100 cursor-pointer box-border b-2 transition-210"
      :class="[dataLabels.enabled 
      ? 'bg-indigo-400 b-indigo-300 hover:(bg-indigo-300 b-indigo-400)'
      : 'bg-indigo-300 b-indigo-400 hover:(bg-indigo-400 b-indigo-300)']"
      @click="dataLabels.enabled=!dataLabels.enabled">
        <div class="text-xl mb-0.5"
        :class="[dataLabels.enabled ? 'i-logos-atomic-icon' : 'i-logos-atomic-icon?mask text-zinc-100']"></div>
        <span class="text-center text-zinc-100">类型</span>
      </div>
    </div>
    <div class="sortedby"></div>
  </div>
  <grid-page v-if="mols4view" :molList="mols4view" :cols='8' :rows="6" class='h-85vh'/>
</div>
</template>
<style scoped>
.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.21s cubic-bezier(0.4, 0, 0.2, 1);
}

.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}
</style>