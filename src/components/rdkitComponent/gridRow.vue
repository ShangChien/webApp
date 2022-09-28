<script setup lang="ts">
import { refDebounced,useElementSize,useVirtualList } from '@vueuse/core'
import { computed, ref,onMounted  } from "vue";
import type {Ref} from "vue";
import type { molData } from "@/components/types";
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import { NInputNumber } from "naive-ui";
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
const input=ref(6)
const cols = refDebounced(input, 1000)
const initArray:molData[] = Array.from(Array(200),(v,i)=>{
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


const postArray:any=computed(()=>{
  return initArray.map((v,i,a)=>{
           let seed=i%cols.value
           if (seed===0){
             return a.slice(i,i+cols.value)
           }
         }).filter((v)=>v!==undefined)
})

const rows=computed(()=>postArray.value.length)

const index: Ref = ref()
const box = ref<HTMLElement> ()
const { width }=useElementSize(box)
const { list, containerProps, wrapperProps, scrollTo }:any = useVirtualList(
  postArray,
  {
    itemHeight: () => (width.value/cols.value + 100),
    overscan: 2,
  },
)
const handleScrollTo = () => {
  scrollTo(index.value)
}
onMounted(() => {
  setTimeout(() => {
      console.log("postlist",postArray.value)
      console.log("outlist",list.value)
  }, 1000);

})
</script>

<template>
<div>
  <n-input-number v-model:value="input" 
                  :update-value-on-input="false"
                  :min="2"
                  class="mb-2 w-20 " 
                  button-placement="both" />
  <div class="inline-block mr-4">
    Jump to index
    <input v-model="index" placeholder="Index" type="number">
  </div>
  <button type="button" @click="handleScrollTo">
    Go
  </button>
  <div ref="box" 
       v-bind="containerProps"
       class="h-78vh overflow-auto p-2 bg-gray-500/5 rounded">
    <div v-bind="wrapperProps">
      <div class="wrapper2" 
        v-for="(data,index) in list" 
        :key="index" 
      >
        <card-rdkit class="w-100\% h-100\%"
          v-for="(itemInner,indexInner) in data.data"
          v-bind="itemInner" 
          :key="indexInner"
        />
      </div>
    </div>
  </div> 
</div>
</template>
<style scoped>
.wrapper2 {
  display: grid;
  grid-template-columns: repeat( v-bind('cols'),1fr);
  grid-column-gap: 0.5em;
}
</style>
