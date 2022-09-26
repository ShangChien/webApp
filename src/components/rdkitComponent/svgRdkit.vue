<script setup lang="ts">
import { ref, onMounted, onUnmounted,watch } from "vue";
import type { molData } from "@/components/types";
const props = defineProps<molData>();
const myWorker = new SharedWorker(new URL('../../worker/sharedWorker.js',import.meta.url),{
  name: 'vastLabSharedWorker',
  type: "module",
})
const bgStyle= ref()
myWorker.port.onmessage = async (e:any)=>{
  requestAnimationFrame(async()=>{
    //console.log(e.data)
    bgStyle.value={
      background: e.data,
      backgroundColor: 'transparent',
      backgroundSize: ['100%','100%'],
    }
    //svg.value.innerHTML = svgText.value// URL.createObjectURL(new Blob([e.data], {type:'image/svg+xml'}))
    //console.log(svgText.value)
  });  
}

myWorker.port.onmessageerror = (e:any)=>{
  console.log('message error:',e)
}
myWorker.onerror = (e:any)=>{
  console.log('error:',e)
}

onMounted(() => {
  //console.log(props)
  myWorker.port.postMessage(JSON.stringify(props))
})

watch(props, (newVal) => {
  myWorker.port.postMessage(JSON.stringify(newVal))
})

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
</style>
