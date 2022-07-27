<script setup lang="ts">
import { ref, onMounted, onUnmounted,watch } from "vue";
import type { molData } from "@/components/types";

const props = defineProps<molData>();
const svg = ref()

let myWorker = new SharedWorker('/src/worker/sharedWorker.js')

//const myWorker:any = inject('myWorker')
myWorker.port.onmessage = (e:any)=>{
  //await nextTick()
  svg.value.innerHTML = e.data; 
}
myWorker.port.onmessageerror = (e:any)=>{
  console.log(e)
}
myWorker.onerror = (e:any)=>{
  console.log(e)
}
onMounted(() => {
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
  <div v-once ><svg ref='svg' viewBox="0 0 200 200" ></svg></div>
  <!-- <canvas ref="canvas" :width="props.width ?? 200" :height="props.height ?? 200" ></canvas> -->
  <!-- <svg style="position:relative;width:100%" v-bind="svgitem.svg">
    <rect v-bind="svgitem.rect" />
    <path
      v-for="item in svgitem.path.hightBonds"
      v-bind="item.path"
    />
    <ellipse
      v-for="item in svgitem.ellipse"
      v-bind="item.ellipse"
    />
    <path
      v-for="item in svgitem.path.symble"
      v-bind="item.path"
    />
    <path
      v-for="item in svgitem.path.bond"
      v-bind="item.path"
    />
  </svg> -->
</template>
<style>
</style>
