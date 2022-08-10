<script setup lang="ts">
import { ref, onMounted, onUnmounted,watch } from "vue";
import type { molData } from "@/components/types";
import { useCopyPNG } from "@/components/rdkitComponent/useCopyPNG";

const myWorker = new SharedWorker(new URL('../../worker/sharedWorker.js',import.meta.url))
const props = defineProps<molData>();
const svg = ref()
const svgText = ref()


myWorker.port.onmessage = async (e:any)=>{
  //await nextTick()
  svg.value.innerHTML = e.data
  svgText.value=e.data
}

myWorker.port.onmessageerror = (e:any)=>{
  console.log(e)
}
myWorker.onerror = (e:any)=>{
  console.log(e)
}

onMounted(() => {
  console.log(myWorker)
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
  <svg ref='svg' viewBox='0, 0, 200, 200' ></svg>
  <button @click="useCopyPNG(svgText)"></button>
</div>
</template>
<style>
</style>
