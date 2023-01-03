<script setup lang="ts">
import { ref, onMounted, onUnmounted,watch,toRaw } from "vue";
import type { molData } from "@/components/types";
const props = defineProps<molData>();

const myWorker = new SharedWorker(new URL('../../worker/sharedWorker.js',import.meta.url),{
  name: 'vastLabSharedWorker',
  type: "module",
})

const bgStyle: any = ref({
      background:null,
      width:"100%",
      height: "100%",
      backgroundColor: 'transparent',
      backgroundSize: ["100%","100%"],
    })
myWorker.port.onmessage = async (e:any)=>{
  requestAnimationFrame(()=>{
    //console.log(e.data)
    bgStyle.value['background'] = e.data
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

myWorker.port.postMessage(JSON.stringify(props))

watch(props, (newVal) => {
  bgStyle.value['background'] = null
  myWorker.port.postMessage(JSON.stringify(newVal))
})

onUnmounted(()=>{
  myWorker.port.close()
})

</script>

<template >
<div>
  <svg v-if="!bgStyle.background" 
       :style='{width:"100%", height:"100%"}'
       viewBox="0,0,200,200"
       class="i-eos-icons-atom-electron c-blue-200"
       ></svg >
  <svg v-else viewBox="0,0,200,200" :style="bgStyle"></svg>
</div>
</template>
<style>
</style>
