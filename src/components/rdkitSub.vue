<script setup lang="ts">
import { ref, onMounted, watch } from 'vue';
import type { molData } from '@/components/types'
import  initRDKitModule from "@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js"
const props = defineProps<molData>()
const rdkitdiv=ref<HTMLDivElement|any>()


function renderMol(props:molData){
  initRDKitModule().then((instance:any)=>{
  const RDKit= instance
  let mol= RDKit.get_mol(props.smiles ?? '')
  let qmol=RDKit.get_qmol(props.qsmiles ?? '')
  console.log(props.smiles)
  let mDetail= JSON.parse(mol.get_substruct_match(qmol))
  mDetail['atoms']=mDetail.atoms?.concat(props.atoms ?? []) ?? props.atoms
  mDetail['bonds']=mDetail.bonds?.concat(props.bonds ?? []) ?? props.bonds
  mDetail['addAtomIndices']=props.addAtomIndices
  mDetail['addBondIndices']=props.addBondIndices
  mDetail['lenged']=props.legend
  mDetail['width']=props.width ?? 100
  mDetail['height']=props.height ?? 100
  mDetail['highlightColour']=props.highlightColor ?? [255,20,147].map((i)=>(i/255))
  mDetail['bondLineWidth']=props.bondLineWidth ?? 1
  mDetail['highlightBondWidthMultiplier']=props.highlightBondWidthMultiplier ?? 28
  mDetail['highlightRadius']=props.highlightRadius ?? 0.4
  mDetail['minFontSize']=props.minFontSize ?? 10
  mDetail['explicitMethyl']=props.explicitMethyl ?? false
  mDetail=JSON.stringify(mDetail)
  let svg=mol.get_svg_with_highlights(mDetail)
  rdkitdiv.value.innerHTML=svg
  //节省内存
  mol.delete()
  qmol.delete()
  })
}
onMounted(()=>{
 renderMol(props)
})
watch(
  props,
  (newVal,oldVal)=>{
    renderMol(newVal)
  }
)

</script>

<template>
<div ref="rdkitdiv"></div>
</template>