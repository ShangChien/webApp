<script setup lang="ts">
import type { molData } from "@/components/types";
import gridPage from "@/components/rdkitComponent/gridPage.vue";
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

</script>

<template>
<div>
  <grid-page :molList="initArray" :cols='8'/>
</div>
</template>
<style scoped>

</style>