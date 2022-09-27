<script setup lang="ts">
import { computed, ref, onMounted } from "vue";
import type { molData } from "@/components/types";
import vGridSvg from "@/components/rdkitComponent/vGridSvg.vue";
import { NButton } from "naive-ui";

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
const colorNum = ref<number>(0);

function change() {
  let color1 = [0.94, 0.475, 0.8];
  let color2 = [0.24, 0.675, 0.8];
  colorNum.value++;
  if (colorNum.value % 2 == 0) {
    mol1.highlightColor = color1;
    mol2.highlightColor = color1;
    mol3.highlightColor = color1;
    mol4.highlightColor = color1;
  } else {
    mol1.highlightColor = color2;
    mol2.highlightColor = color2;
    mol3.highlightColor = color2;
    mol4.highlightColor = color2;
  }
}
let col = 6
const initArray = Array.from(Array(100),(v,i)=>{
  let seed=i%6
  if (seed===0){
    return mol1
  }else if (seed===1){
    return mol2
  }else if (seed===2){
    return mol3
  }else if (seed===3){  
    return mol4
  }else if (seed===4){
    return mol5
  }else if (seed===5){
    return mol6
  }
})

const postArray:any=initArray.map((v,i,a)=>{
  let seed=i%col
  if (seed===0){
    return a.slice(i,i+6)
  }
}).filter((v)=>v!==undefined)
onMounted(()=>{
  console.log(postArray)
})
</script>

<template>
<div>
  <n-button @click="change"> 切换颜色 </n-button><br>
  <div class="wrapper1">
    <v-grid-svg
      v-for="(item,index) in postArray"
      :molList="item" 
      :key="index"
    />
  </div>
</div>
</template>
<style scoped>
.wrapper1 {
  display: grid;
  grid-template-columns: 1fr;
  grid-row-gap: 0.5em;
}
</style>
