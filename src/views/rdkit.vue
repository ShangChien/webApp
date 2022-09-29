<script setup lang="ts">
import { refDebounced } from '@vueuse/core'
import { computed, ref, onMounted,onUpdated,onBeforeMount } from "vue";
import type { molData } from "@/components/types";
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import { NInputNumber,NPagination,NPopover,NTag } from "naive-ui";
const mol1: molData = {
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  qsmiles: "CC(=O)Oc1ccccc1C(=O)O",  
};
const mol2: molData = {
  smiles:
    "CSCC[C@H](NC(=O)[C@H](CC1=CNC2=C1C=CC=C2)NC(=O)CCNC(=O)OC(C)(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC1=CC=CC=C1)C(N)=O",
  qsmiles: "CS.c1ccccc1.c1ccccc1",
  //highlightColor: [0.624, 0.675, 0.2],
};
const mol3: molData = {
  smiles:
    "CC(C)(C)C(C=C1)=CC2=C1SC3=C2N(C4=CC(C5=CC=C(C(C)(C)C)C=C5)=CC6=C4B3C7=C8N6C9=C(C=CC=C9)C8=CC%10=C7C=CC=C%10)C%11=CC=C(C(C)(C)C)C=C%11",
  qsmiles: "[B].c1ccccc1",
  atoms:{3:[2,4,5,8],5:[12,32,23]},
  bonds:{3:[2,4,5,8]}
};
const mol4: molData = {
  smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)",
  qsmiles: "[O].[O]"
};
const mol5: molData = {
  smiles:
    "CC(C)(C)C(C=C1)=CC2=C1SC3=C2N(C4=CC(C5=CC=C(C(C)(C)C)C=C5)=CC6=C4B3C7=C8N6C9=C(C=CC=C9)C8=CC%10=C7C=CC=C%10)C%11=CC=C(C(C)(C)C)C=C%11",
  qsmiles: "[B].c1ccccc1",
  atoms:{3:[2,4,5,8],5:[12,32,23]},
  bonds:{3:[2,4,5,8]}
};
const mol6: molData = {
  smiles:
    "CC(C)(C)C(C=C1)=CC2=C1SC3=C2N(C4=CC(C5=CC=C(C(C)(C)C)C=C5)=CC6=C4B3C7=C8N6C9=C(C=CC=C9)C8=CC%10=C7C=CC=C%10)C%11=CC=C(C(C)(C)C)C=C%11",
  qsmiles: "[B].c1ccccc1",
  atoms:{3:[2,4,5,8],5:[12,32,23]},
  bonds:{3:[2,4,5,8]}
};

const initArray:molData[] = Array.from(Array(20000),(v,i)=>{
  let seed=i%7
  if (seed===1){
    return mol1
  }else if (seed===2){
    return mol2
  }else if (seed===3){
    return mol3
  }else if (seed===4){  
    return mol4
  }else if (seed===5){
    return mol5
  }else{
    return mol6
  }
})

const inputCols = ref(8)
const inputRows = ref(8)
const cols = refDebounced(inputCols, 500)
const rows = refDebounced(inputRows, 500)
const pageSize = computed(()=>cols.value*rows.value)
const pageCount = computed(()=>Math.ceil(initArray.length/pageSize.value))
const currentPageInput=ref(1)
const currentPage = refDebounced(currentPageInput, 500)
const arrayShow=computed(()=>{
  let start = (currentPage.value-1)*pageSize.value
  let end = currentPage.value*pageSize.value
  return initArray.slice(start,end)
})
// const postArray:any=computed(()=>{
//   return initArray.map((v,i,a)=>{
//            let seed=i%cols.value
//            if (seed===0){
//              return a.slice(i,i+cols.value)
//            }
//          }).filter((v)=>v!==undefined)
// })

</script>

<template>
<div>
  <div class="grid grid-cols-2">
    <div>
      <n-pagination v-model:page="currentPageInput"
                    :page-count="pageCount"
                    :page-slot="7"
                    class="mt--2 mb-2 "
                    size="medium"
                    show-quick-jumper>
        <template #goto>
          跳至:
        </template>
      </n-pagination>
    </div>
    <div class="justify-self-end mt--2 mb-2">
      <n-tag :bordered="false" class="mr-1 bg-red-50">
        <span class="text-1.2em ">Total: {{initArray.length}}</span>
      </n-tag>
      <n-tag :bordered="false" class="mr-1 bg-red-50">
        <span class="text-1.2em "> ( {{pageSize}} / page)</span>
      </n-tag>
      <n-popover trigger="click">
        <template #trigger>
          <div class="i-fluent-table-settings-20-filled text-2.5em c-indigo-400 mt--1.2"></div>
        </template>
        <template #default>
          <div class="grid justify-items-end grid-cols-2 gap-2">
            <span class=" pr-1 pl-2 text-xl">cols:</span>
            <n-input-number v-model:value="inputCols" 
                            :update-value-on-input="false"
                            :min="2"
                            size="small"
                            class="w-20 " 
                            button-placement="both" />
            <span class=" pr-1 pl-2 text-xl ">rows:</span>
            <n-input-number v-model:value="inputRows" 
                            :update-value-on-input="false"
                            :min="2"
                            size="small"
                            class="w-20" 
                            button-placement="both" />
          </div>
        </template>
      </n-popover>
    </div>
  </div>
  <div class="wrapper1" >
      <card-rdkit class="w-100\% h-100\%"
        v-for="(itemInner,indexInner) of arrayShow"
        v-bind="itemInner" 
        :key="indexInner"
      />
  </div>
</div>
</template>
<style scoped>
.wrapper1 {
  display: grid;
  grid-template-columns: repeat( v-bind('cols'),1fr);
  grid-template-rows: repeat(v-bind(rows),1fr);
  grid-row-gap: 0.5em;
  grid-column-gap: 0.5em;
}
</style>