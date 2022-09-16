<script setup lang="ts">
import { ref, onMounted, onUnmounted,watch,computed,reactive,nextTick } from "vue";
import type { molData } from "@/components/types";
import { optimize } from 'svgo/lib/svgo.js';
const myWorker = new SharedWorker(new URL('../../worker/sharedWorker.js',import.meta.url))
const props = defineProps<molData>();
const res = ref(false)
const svgText = ref()
function cssBgSvg(svgI: string) {
  return 'data:image/svg+xml;utf8,'+ optimize(svgI, {}).data
    .replace('<svg', (~svgI.indexOf('xmlns') ? '<svg' : '<svg xmlns="http://www.w3.org/2000/svg"'))
    .replace(/"/g, '\'')
    .replace(/%/g, '%25')
    .replace(/#/g, '%23')
    .replace(/{/g, '%7B')
    .replace(/}/g, '%7D')
    .replace(/</g, '%3C')
    .replace(/>/g, '%3E')
    .replace(/atom/g,'a')
    .replace(/bond/g,'b')
}
const bgSvg = computed(() => res.value ? 'url("'+cssBgSvg(svgText.value)+'")': 'url("")')
const bgStyle= ref()
//const bgStyle:any=computed(()=>{ Res.value ? `background:url("${cssBgSvg(svgText)}") no-repeat center;` : '' })

myWorker.port.onmessage = async (e:any)=>{
  svgText.value = e.data
  //requestAnimationFrame(()=>{
    //svg.value.innerHTML = svgText.value// URL.createObjectURL(new Blob([e.data], {type:'image/svg+xml'}))
    //console.log(svgText.value)
  //}); 
  res.value=true
  bgStyle.value={
    background: bgSvg,
    backgroundColor: 'transparent',
    backgroundSize: ['100%','100%'],
  }
}

myWorker.port.onmessageerror = (e:any)=>{
  console.log(e)
}
myWorker.onerror = (e:any)=>{
  console.log(e)
}

onMounted(() => {
  //console.log(props)
  myWorker.port.postMessage(JSON.stringify(props))
})

// watch(props, (newVal) => {
//   myWorker.port.postMessage(JSON.stringify(newVal))
// })

onUnmounted(()=>{
  myWorker.port.close()
})
</script>

<template >
<div>
  <svg viewBox="0,0,200,200" :style="bgStyle"></svg>
</div>
</template>
<style>
.svgMol {
  background-color: transparent;
  background-size: 100% 100%;
}
</style>
