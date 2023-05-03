<script setup lang="ts">
import { computed,ref,reactive } from 'vue'
import { NButton } from 'naive-ui'
import { createReusableTemplate } from '@vueuse/core'
import { useMolStore } from '@/stores/molStore'
import type { molData } from "@/components/types"
import gridPage from "@/components/rdkitComponent/gridPage.vue";
interface dataFilter {
  items: string[],
  selected: string[],
  expend:boolean
}
const [DefineTag, ReuseTag] = createReusableTemplate<{label:string }>()
const store = useMolStore()
const labels = computed(()=>store.getAllLabels)
const dataLabels:dataFilter=reactive({items:labels,selected:[],expend:false})
const dataTypes:dataFilter=reactive({items:[],selected:[],expend:false})
const mols4view = computed(()=>{
  return store.getByTypesAndLabels(dataTypes.selected,dataLabels.selected)
})
function btclick(label:string) {
  return dataLabels.selected.includes(label) ? dataLabels.selected.splice(dataLabels.selected.indexOf(label), 1) : dataLabels.selected.push(label)
}

</script>

<template>
<div>
  <DefineTag v-slot="{label}">
    <n-button class="m-1" strong secondary size='small'
      :type="dataLabels.selected.includes(label) ? 'primary':'tertiary'"
      @click="btclick(label)">
      {{label}}
    </n-button>
  </DefineTag>
  <div class="filter">
    <div class="flex flex-nowrap justify-between items-center">
      <div class="flex justify-begin overflow-hidden flex-auto "
        :class="[dataLabels.expend ? 'flex-nowrap':'flex-wrap']">
        <ReuseTag v-for="label in labels" :label="label"></ReuseTag>
      </div>
      <div class="rd-2 m-2 bg-blue-100">
        <div @click="dataLabels.expend = !dataLabels.expend"
        :class="[dataLabels.expend ? 'i-ion-chevron-expand' : 'i-ion-chevron-collapse']"
        class="flex-none  text-2xl  c-zinc-400 hover:c-slate-600" />
      </div>
    </div>
    <div class="type"></div>
    <div class="sortedby"></div>
  </div>
  <grid-page v-if="mols4view" :molList="mols4view" :cols='8' :rows="6" class='h-85vh'/>
</div>
</template>
<style scoped>

</style>