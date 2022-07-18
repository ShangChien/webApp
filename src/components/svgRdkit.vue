<script setup lang="ts">
import { ref, onMounted, onUnmounted,inject,nextTick } from "vue";
import type { molData } from "@/components/types";

const props = defineProps<molData>();
let timeFn = ref()
const svg = ref()

let myWorker = new SharedWorker('/src/worker/sharedWorker.js')

//const myWorker:any = inject('myWorker')
myWorker.port.onmessage= async (e:any)=>{
  //await nextTick()
  svg.value.innerHTML = e.data; 
  console.log(props.smiles,myWorker)
  myWorker.port.close()
}
myWorker.port.onmessageerror = (e:any)=>{
  console.log(e)
}
myWorker.onerror = (e:any)=>{
  console.log(e)
}
onMounted(() => {
  console.log(props.smiles,myWorker)
  myWorker.port.postMessage(JSON.stringify(props))
})

// watch(props, (newVal) => {
//   let dataNew = toRaw(newVal)
//   svg.value.innerHTML=renderMol(dataNew);
// });
// watch(
//   data,
//   async (dataNew) => {
//     terminate()
//     await nextTick(()=>{
//       requestAnimationFrame(() => {
//         svg.value.innerHTML=dataNew;
//       })
//     })
    
//   }
// )

onUnmounted(()=>{
  clearTimeout(timeFn.value)
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
