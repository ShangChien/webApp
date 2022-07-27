<script setup lang="ts">
import { ref, onMounted, watch } from "vue";
import type { molData } from "@/components/types";
import initRDKit from "@/components/rdkitComponent/RDKit";
const props = defineProps<molData>();
const rdkitdiv = ref<HTMLDivElement | any>();

function renderMol(props: molData) {
  window.RDKit.prefer_coordgen(true);
  let mol = window.RDKit.get_mol(props.smiles ?? "");
  let qmol = window.RDKit.get_qmol(props.qsmiles ?? "");
  let mDetail = JSON.parse(mol.get_substruct_match(qmol));
  mDetail["atoms"] = mDetail.atoms?.concat(props.atoms ?? []) ?? props.atoms;
  mDetail["bonds"] = mDetail.bonds?.concat(props.bonds ?? []) ?? props.bonds;
  mDetail["addAtomIndices"] = props.addAtomIndices;
  mDetail["addBondIndices"] = props.addBondIndices;
  mDetail["lenged"] = props.legend;
  mDetail["width"] = props.width ?? 200;
  mDetail["height"] = props.height ?? 200;
  mDetail["highlightColour"] =
    props.highlightColor ?? [208, 78, 214].map((i) => i / 255);
  mDetail["bondLineWidth"] = props.bondLineWidth ?? 1;
  mDetail["highlightBondWidthMultiplier"] =
    props.highlightBondWidthMultiplier ?? 28;
  mDetail["highlightRadius"] = props.highlightRadius ?? 0.4;
  mDetail["minFontSize"] = props.minFontSize ?? 10;
  mDetail["explicitMethyl"] = props.explicitMethyl ?? false;
  mDetail = JSON.stringify(mDetail);
  let svg = mol.get_svg_with_highlights(mDetail);
  //console.log(svg)
  rdkitdiv.value.innerHTML = svg;
  //节省内存
  mol.delete();
  qmol.delete();
}
onMounted(async () => {
  await initRDKit.then((res) => {
    window.RDKit = res;
  });
  renderMol(props);
});
watch(props, (newVal) => {
  renderMol(newVal);
});
</script>
<template>
  <div ref="rdkitdiv" class="svg"></div>
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
