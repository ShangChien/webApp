<script setup lang="ts">
import { ref, watch,computed, reactive,inject,toRaw } from "vue";
import type { Ref } from 'vue'
import { useRefHistory } from '@vueuse/core'
import type { molData } from "@/components/types";
import { useGetSvg } from '@/components/rdkitComponent/composable/useGetSvg'
const props = defineProps<molData>();
const rdkit = inject('rdkit')
const svgItem = ref(useGetSvg(props,rdkit))
const highlightMap:Ref<molData>|any = ref({
  id: 0,
  smiles: props.smiles,
  atoms: props.atoms,
  bonds: props.bonds,
  label: props.labels,
});
const { history, undo, redo, clear,canUndo, canRedo } = useRefHistory(highlightMap, { deep: true })
defineExpose({ history, undo, redo, clear, canUndo, canRedo, highlightMap })
const siteType:any= inject('siteType')
const Color = (n:any) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,70%,1)'
//目的,遍历highlightMap对象atoms属性的所有子属性，如果包含atomIndex，则删除
function preHandleIndex(obj:object|undefined,atomIndex:number){
  highlightMap.value.id++
  for (let key in obj) {
    if (obj[key].includes(atomIndex)) {
      obj[key].splice(obj[key].indexOf(atomIndex),1);
      if (obj[key].length === 0){
        delete obj[key]
      }
    }
  }
}
function domClick($event: any) {
  console.log('before',props)
  let itemList = $event.target.getAttribute("class").split(" ")[0].split("-");
  if (itemList[0] == "atom") {
    //如果点击atom
    //改颜色
    //$event.target.style.stroke = Color(siteType.value);
    //$event.target.style.fill = Color(siteType.value);
    //$event.target.style.opacity = 0.7;
    var atomIndex = itemList[1] / 1
    preHandleIndex(highlightMap.value.atoms, atomIndex)
    //添加index到数组
    highlightMap.value.atoms[siteType.value] = highlightMap.value.atoms[siteType.value] ? 
                                         highlightMap.value.atoms[siteType.value]:[]
    highlightMap.value.atoms[siteType.value].push(itemList[1] / 1);
    //highlightMap.atoms[siteType.value] = Array.from(new Set(highlightMap.atoms[siteType.value])).sort();
    //emit("update-mol", highlightMap);
    console.log('after',props)
  } else if (itemList[0] == "bond") {
    //如果点击bond
    //$event.target.style.stroke = Color(siteType.value);
    //$event.target.style.fill = Color(siteType.value);
    //$event.target.style.opacity = 0.6;
    var bondIndex = itemList[1] / 1
    preHandleIndex(highlightMap.value.bonds, bondIndex)
    //添加index到数组
    highlightMap.value.bonds[siteType.value] = highlightMap.value.bonds[siteType.value] ? 
                                         highlightMap.value.bonds[siteType.value]:[]
    highlightMap.value.bonds[siteType.value].push(itemList[1] / 1);
    //highlightMap.bonds[siteType.value] = Array.from(new Set(highlightMap.bonds[siteType.value])).sort();
    //emit("update-mol", highlightMap);
    //console.log(highlightMap.bonds)
  } else {
    console.log(itemList[0], "error");
  }
}
function domDblClick($event: any) {
  let itemList = $event.target.getAttribute("class").split(" ")[0].split("-");
  if (itemList[0] == "atom") {
    //如果点击atom
    //改颜色
    //$event.target.style.stroke = "#9FACE6";
    //$event.target.style.fill = "#9FACE6";
    //$event.target.style.opacity = 0.3;
    //删除index
    preHandleIndex(highlightMap.value.atoms, itemList[1]/1)
    //emit("update-mol", highlightMap);
    //console.log(highlightMap.highlightAtoms)
  } else if (itemList[0] == "bond") {
    //如果点击bond
    //$event.target.style.stroke = "#9FACE6";
    //$event.target.style.fill = "#9FACE6";
    //$event.target.style.opacity = 0.3;
    //删除index
    preHandleIndex(highlightMap.value.bonds, itemList[1]/1)
    //emit("update-mol", highlightMap);
    //console.log(highlightMap.highlightAtoms)
  } else {
    console.log(itemList[0], "error");
  }
}
function clearAll() {
  highlightMap.value.atoms=highlightMap.value.atoms.forEach()
  highlightMap.value.bonds={}
  svgItem.value = useGetSvg(props,rdkit)
  highlightMap.id=0
  //emit("update-mol", highlightMap);
}

//根据highlightMap初始化高亮svg中的atoms和bonds
//const svg_id = ref()
// function initHighlightSvg() {
//   //初始化处理高亮atoms
//   for (let item of svg_id.value.getElementsByTagName("ellipse")) {
//     //在atoms中添加item.class属性具有highlightMap.atoms的元素
//     if (
//       highlightMap.atoms?.includes(item.getAttribute("class").split("-")[1] / 1)
//     ) {
//       item.style.stroke = Color(siteType.value);
//       item.style.fill = Color(siteType.value);
//       item.style.opacity = 0.8;
//     }
//   }
//   //初始化处理高亮bonds
//   for (let item of svg_id.value.getElementsByTagName("path")) {
//     if (
//       item.style.strokeWidth === "13.7px" &&
//       highlightMap.bonds?.includes(
//         item.getAttribute("class").split(" ")[0].split("-")[1] / 1
//       )
//     ) {
//       item.style.stroke = Color(siteType.value);
//       item.style.fill = Color(siteType.value);
//       item.style.opacity = 0.4;
//     }
//   }
// }

watch(props, (newVal) => {
  svgItem.value = useGetSvg(newVal,rdkit)
  highlightMap.value.smiles = newVal.smiles;
  highlightMap.value.atoms = newVal.atoms;
  highlightMap.value.bonds = newVal.bonds;
  console.log("watch", highlightMap.value);
});
</script>

<template class="svg">
  <svg :key="highlightMap.id"
    v-bind="svgItem.svg"
    class="svgstyle"
    @click.right.prevent="clearAll"
    style="width: 100%; height: 100%"
  >
    <rect v-bind="svgItem.rect" />
    <path
      v-for="item in svgItem.path.symble"
      v-bind="item.path"
    />
    <path
      v-for="item in svgItem.path.bond"
      v-bind="item.path"
    />
    <path
      v-for="item in svgItem.path.hightBonds"
      v-bind="item.path"
      @click="domClick($event)"
      @dblclick="domDblClick($event)"
    />
    <ellipse
      v-for="item in svgItem.ellipse"
      v-bind="item.ellipse"
      @click="domClick($event)"
      @dblclick="domDblClick($event)"
    />
  </svg>
</template>
<style>
.svg {
  display: inline-block;
  position: relative;
  width: 100%;
  padding-bottom: 100%;
  vertical-align: middle;
  overflow: hidden;
}
</style>
