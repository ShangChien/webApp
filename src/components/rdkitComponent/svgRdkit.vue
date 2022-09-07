<script setup lang="ts">
import { ref, onMounted, onUnmounted,watch,toRaw } from "vue";
import type { molData } from "@/components/types";

const myWorker = new SharedWorker(new URL('../../worker/sharedWorker.js',import.meta.url))
const props = defineProps<molData>();
const svg = ref()
const svgText = ref()
function preHandleProps(props:molData) {
  const concatIndex = (list1:{[key: number|string]: number[]}|undefined) =>{ 
    let outList: any[] = []
    for(let i in list1){
      outList = outList.concat(list1[i])
    }
    return Array.from(new Set(outList)) 
  }
  const initProps:any = toRaw(props)
  initProps.atoms = concatIndex(props.atoms)
  initProps.bonds = concatIndex(props.bonds)
  return initProps
}




myWorker.port.onmessage = async (e:any)=>{
  //await nextTick()
  requestAnimationFrame(()=>{
    svg.value.innerHTML = e.data// URL.createObjectURL(new Blob([e.data], {type:'image/svg+xml'}))
  }); 
  svgText.value=e.data
}

myWorker.port.onmessageerror = (e:any)=>{
  console.log(e)
}
myWorker.onerror = (e:any)=>{
  console.log(e)
}

onMounted(() => {
  //console.log(props)
  myWorker.port.postMessage(JSON.stringify(preHandleProps(props)))
})

watch(props, (newVal) => {
  myWorker.port.postMessage(JSON.stringify(preHandleProps(newVal)))
})

onUnmounted(()=>{
  myWorker.port.close()
})
</script>

<template >
<div>
  <svg ref='svg' viewBox='0, 0, 200, 200' ></svg>
</div>
</template>
<style>
</style>
