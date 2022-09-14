<script setup lang="ts">
import { ref, watch, onMounted, computed,nextTick, inject } from "vue";
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
const highlightMap:Ref< molData>|any =ref({
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
  console.log(canRedo.value,canUndo.value)
  if (canUndo.value ) {
    undo()
    svgItem.value=useGetSvg(highlightMap.value,rdkit)
  } 
}
function redoRender(){
  console.log(canRedo.value,canUndo.value)
  if (canRedo.value ) {
    redo()
    svgItem.value=useGetSvg(highlightMap.value,rdkit)
  } 
}
defineExpose({ history, undoRender, redoRender, clear, canUndo, canRedo, highlightMap })
//侦听props，截取props初始值，刷新渲染
watch(
  props,
  async (newVal,oldVal)=>{
    let { qsmiles,smiles,atoms,bonds,labels } = JSON.parse(JSON.stringify(newVal))
    highlightMap.value = {
      id: 0,
      qsmiles: qsmiles,
      smiles: smiles,
      atoms: atoms,
      bonds: bonds,
      label: labels,
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
let timer:any
function domClick($event: any) {
  
  clearTimeout(timer)
  timer = setTimeout(() => {
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
  }, 200);
}
function domDblClick($event: any) {
  clearTimeout(timer);
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
function clearAll() {
  //await nextTick()
  console.log('clearAll',history.value,initState.value)
  highlightMap.value = JSON.parse(JSON.stringify(initState.value))
  svgItem.value = useGetSvg(initState.value, rdkit)
  //emit("update-mol", highlightMap);
}
onMounted(()=>{
console.log(canRedo.value,canUndo.value)
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

<template class="svg">
  <svg 
    v-bind="svgItem.svg"
    class="svgstyle"
    @click.right.prevent="clearAll"
    style="width: 100%; height: 100%"
  >
    <rect v-bind="svgItem.rect" />
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
      @dblclick="domDblClick($event)"
    />
    <ellipse
      v-for="(item,index) in svgItem.ellipse"
      v-bind="item.ellipse"
      :key="index"
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
