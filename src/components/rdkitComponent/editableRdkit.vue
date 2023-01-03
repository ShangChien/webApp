<script setup lang="ts">
import { ref, watch, onMounted, computed,onBeforeUpdate, nextTick, inject } from "vue";
import type { Ref } from 'vue'
import { useRefHistory } from '@vueuse/core'
import type { molData } from "@/components/types";
import { useGetSvg } from '@/components/rdkitComponent/composable/useGetSvg'
//inject global state
const rdkit = inject('rdkit')
const siteType:any= inject('siteType')
const Color = (n:any) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,70%,1)'
//截取初始的props
const props = defineProps<molData>();
let { qsmiles,smiles,atoms,bonds,labels } = JSON.parse(JSON.stringify(props)) 
const highlightMap:Ref< molData> =ref({
  id: 0,
  qsmiles: qsmiles,
  smiles: smiles,
  atoms: atoms,
  bonds: bonds,
  label: labels,
});
const initState=ref(JSON.parse(JSON.stringify(highlightMap.value)))
const svgItem = ref(useGetSvg(highlightMap.value,rdkit))
//记录历史状态
const { history, undo, redo, clear, canUndo, canRedo } = useRefHistory(highlightMap, { deep: true, flush: 'post' })

function undoRender(){
  if (canUndo.value) {
    undo()
    svgItem.value=useGetSvg(highlightMap.value,rdkit)
  } 
}
function redoRender(){
  if (canRedo.value ) {
    redo()
    svgItem.value=useGetSvg(highlightMap.value,rdkit)
  } 
}
defineExpose({ history, undoRender, redoRender, canUndo, canRedo, highlightMap })
//侦听props，截取props初始值，刷新渲染
// onBeforeUpdate(()=>{
//   svgItem.value = useGetSvg(highlightMap.value,rdkit)
// })
watch(
  props,
  async (newVal)=>{
    let { qsmiles,smiles,atoms,bonds,labels } = JSON.parse(JSON.stringify(newVal))
    highlightMap.value = {
      id: 0,
      qsmiles: qsmiles,
      smiles: smiles,
      atoms: atoms,
      bonds: bonds,
      labels: labels,
    }
    await nextTick()
    clear()
    initState.value = JSON.parse(JSON.stringify(highlightMap.value))
    svgItem.value = useGetSvg(highlightMap.value,rdkit)
  }
)

//目的,遍历highlightMap对象atoms属性的所有子属性，如果包含atomIndex，则删除
function preHandleIndex(obj:object|undefined, atomIndex:number){
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
  
  let itemList = $event.target.getAttribute("class").split(" ")[0].split("-");
  if (itemList[0] == "atom") {
    //如果点击atom
    //改颜色
    $event.target.style.stroke = Color(siteType.value);
    $event.target.style.fill = Color(siteType.value);
    $event.target.style.opacity = 0.7;
    var atomIndex = itemList[1] / 1
    preHandleIndex(highlightMap.value.atoms, atomIndex)
    //添加index到数组
    highlightMap.value.atoms[siteType.value] = highlightMap.value.atoms[siteType.value] ? 
                                         highlightMap.value.atoms[siteType.value]:[]
    highlightMap.value.atoms[siteType.value].push(itemList[1] / 1);
    //highlightMap.atoms[siteType.value] = Array.from(new Set(highlightMap.atoms[siteType.value])).sort();
    //emit("update-mol", highlightMap);
  } else if (itemList[0] == "bond") {
    //如果点击bond
    $event.target.style.stroke = Color(siteType.value);
    $event.target.style.fill = Color(siteType.value);
    $event.target.style.opacity = 0.6;
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

function domRClick($event: any) {
  let itemList = $event.target.getAttribute("class").split(" ")[0].split("-");
  if (itemList[0] == "atom") {
    //如果点击atom
    //改颜色
    $event.target.style.stroke = "#9FACE6";
    $event.target.style.fill = "#9FACE6";
    $event.target.style.opacity = 0.3;
    //删除index
    preHandleIndex(highlightMap.value.atoms, itemList[1]/1)
    //emit("update-mol", highlightMap);
    //console.log(highlightMap.highlightAtoms)
  } else if (itemList[0] == "bond") {
    //如果点击bond
    $event.target.style.stroke = "#9FACE6";
    $event.target.style.fill = "#9FACE6";
    $event.target.style.opacity = 0.3;
    //删除index
    preHandleIndex(highlightMap.value.bonds, itemList[1]/1)
    //emit("update-mol", highlightMap);
    //console.log(highlightMap.highlightAtoms)
  } else {
    console.log(itemList[0], "error");
  }
}

onMounted(()=>{
//console.log(svgItem.value)
})
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
</script>

<template >
  <svg 
    v-bind="svgItem.svg"
    class="w-100% h-100% svg bg-white "
    viewBox="0 0 200 200"
  >
    <path 
      v-for="(item,index) in svgItem.path.symble"
      v-bind="item.path"
      :key="index"
    />
    <path
      v-for="(item,index) in svgItem.path.bond"
      v-bind="item.path"
      :key="index"
    />
    <path
      v-for="(item,index) in svgItem.path.hightBonds"
      v-bind="item.path"
      :key="index"
      @click="domClick($event)"
      @click.right.prevent="domRClick($event)"
    />
    <ellipse 
      v-for="(item,index) in svgItem.ellipse"
      v-bind="item.ellipse"
      :key="index"
      @click="domClick($event)"
      @click.right.prevent="domRClick($event)"
    />
  </svg>
</template>
<style>
.svg {
  display:flex;
  align-items:center;
  justify-content:center;
}
</style>
