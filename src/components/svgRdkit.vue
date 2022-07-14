<script setup lang="ts">
import { ref, onMounted, onUnmounted, nextTick, toRaw,
  inject, watch, reactive } from "vue";
import type { molData } from "@/components/types";
import { useWebWorker } from '@vueuse/core'
//import initRDKit from "@/components/RDKit";
//const { data, post, terminate } = useWebWorker('/src/worker/workerRdkit.js')

//const canvas=ref()
const props = defineProps<molData>();
//const rdkit:any = toRaw(inject('rdkit'))
const svg=ref()

// function renderMol(props:any){
//   rdkit.prefer_coordgen(true);
//   let mol = rdkit.get_mol(props.smiles ?? "");
//   let qmol = rdkit.get_qmol(props.qsmiles ?? "");
//   //let mDetail = JSON.parse(mol.get_substruct_match(qmol))
//   let mdetailsRaw = mol.get_substruct_matches(qmol);
//   let mDetail = mdetailsRaw.length > 2 ? JSON.parse(mdetailsRaw) : [];
//   mDetail = mDetail.reduce(
//     (acc: any, { atoms, bonds }: any) => ({
//       atoms: [...acc.atoms, ...atoms],
//       bonds: [...acc.bonds, ...bonds],
//     }),
//     { bonds: [], atoms: [] }
//   );
//   mDetail["atoms"] = mDetail.atoms?.concat(props.atoms ?? []) ?? props.atoms;
//   mDetail["bonds"] = mDetail.bonds?.concat(props.bonds ?? []) ?? props.bonds;
//   mDetail["addAtomIndices"] = props.addAtomIndices;
//   mDetail["addBondIndices"] = props.addBondIndices;
//   mDetail["lenged"] = props.legend;
//   mDetail["width"] = props.width ?? 200;
//   mDetail["height"] = props.height ?? 200;
//   mDetail["highlightColour"] =
//     props.highlightColor ?? [208, 78, 214].map((i) => i / 255);
//   mDetail["bondLineWidth"] = props.bondLineWidth ?? 1;
//   mDetail["highlightBondWidthMultiplier"] =
//     props.highlightBondWidthMultiplier ?? 15;
//   mDetail["highlightRadius"] = props.highlightRadius ?? 0.25;
//   mDetail["minFontSize"] = props.minFontSize ?? 10;
//   mDetail["explicitMethyl"] = props.explicitMethyl ?? false;
//   mDetail = JSON.stringify(mDetail);
//   let out = mol.get_svg_with_highlights(mDetail);
//   //draw canvas
//   //mol.draw_to_canvas_with_highlights(canvas.value, mDetail);
//   qmol.delete();
//   mol.delete();
//   return out
// }
const myWorker = new SharedWorker('/src/worker/shareWorker.js')
onMounted(async () => {
  const propsData:any = toRaw(props)
  propsData.atoms = toRaw(props.atoms)
  propsData.bonds = toRaw(props.bonds)
  myWorker.port.postMessage(propsData);
  myWorker.port.onmessage = async function(e) {
    await nextTick(()=>{
      requestAnimationFrame(() => {
        svg.value.innerHTML=e.data;
      })
    })  
  }
  //post(propsData)
  //await nextTick()
   // await Promise.resolve().then(()=>{
  //   setTimeout(() => {
  //     renderMol(props)
  //   },0)
  // })
});

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


// onUnmounted(()=>{
//   terminate()
// })
</script>

<template >
  <div ><svg ref='svg' viewBox="0 0 200 200" ></svg></div>
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
