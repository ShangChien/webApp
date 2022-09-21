<script setup lang="ts">
import { reactive, ref } from "vue";
import type { molData } from "@/components/types";
import rdkitSub from "@/components/rdkitComponent/rdkitSub.vue";
import svgRdkit from "@/components/rdkitComponent/svgRdkit.vue";
import { NSpace, NButton,NCollapse,NCollapseItem } from "naive-ui";

const mol1: molData = reactive({
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  qsmiles: "CC(=O)Oc1ccccc1C(=O)O",
  width: 200,
  height: 200,
  addAtomIndices: true,
});

const mol2: molData = reactive({
  smiles:
    "CSCC[C@H](NC(=O)[C@H](CC1=CNC2=C1C=CC=C2)NC(=O)CCNC(=O)OC(C)(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC1=CC=CC=C1)C(N)=O",
  qsmiles: "CS.c1ccccc1.c1ccccc1",
  width: 400,
  height: 400,
  //highlightColor: [0.624, 0.675, 0.2],
  addAtomIndices: true,
  addBondIndices: true,
});

const mol3: molData = reactive({
  smiles:
    "CC(C)(C)C(C=C1)=CC2=C1SC3=C2N(C4=CC(C5=CC=C(C(C)(C)C)C=C5)=CC6=C4B3C7=C8N6C9=C(C=CC=C9)C8=CC%10=C7C=CC=C%10)C%11=CC=C(C(C)(C)C)C=C%11",
  qsmiles: "[B].c1ccccc1",
  atoms:{},
  width: 200,
  height: 200,
  css:true,
  addAtomIndices: false,
  addBondIndices: false,
  highlightColor: [0.24, 0.675, 0.8],
});
const mol4: molData = reactive({
  smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)",
  qsmiles: "[O].[O]",
  width: 200,
  height: 200,
  addBondIndices: true,
  addAtomIndices: true,
  //highlightColor: [0.94, 0.475, 0.8],
});
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
const show=ref<boolean>(true)
</script>

<template>
<div style="position:relative">
  <NButton @click="change"> 切换颜色 </NButton>
  <NButton @click="show=!show"> show </NButton>
   <n-collapse :default-expanded-names="['1']" >
    <n-collapse-item  display-directive="if" name="1" >
    <div style="padding-left:20px;padding-right:20px">
      <svg-rdkit v-bind="mol3"  style="width:20%" ></svg-rdkit>
    </div>
    </n-collapse-item>
  </n-collapse>
  <n-space v-if="show">
  <n-space>
 
    <rdkit-sub v-bind="mol4"></rdkit-sub>
  </n-space>
  <n-space>
    <rdkit-sub v-bind="mol2"></rdkit-sub>
    <rdkit-sub v-bind="mol1"></rdkit-sub>
  </n-space>
  </n-space>
  
</div>
</template>
