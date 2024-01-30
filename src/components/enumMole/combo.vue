<script setup lang='ts'>
import { reactive, ref, watch } from 'vue'
import cardRdkit from '@/components/rdkitComponent/cardRdkit.vue'
import comboSetting from '@/components/enumMole/comboSetting.vue'
import type { molData } from '@/components/types'

// use&dataState
const props = defineProps<{ ligands: molData[] | any; cores: molData[] | any }>()
const coresData = ref()
watch(
  () => props.cores,
  (nCores) => {
    coresData.value = nCores.map((mol: molData) => {
      const Nmol = JSON.parse(JSON.stringify(mol))
      // 必须新建空对象'enumAtoms'，才能深层赋值
      Nmol.enumAtoms = {}
      Object.keys(mol.atoms ?? {}).forEach((colorIndex) => {
        // const step = (100 / mol.atoms[colorIndex].length).toFixed(2)
        Nmol.enumAtoms[colorIndex] = reactive({
          array: mol.atoms[colorIndex],
          range: [0, mol.atoms[colorIndex].length],
          keepSame2Index: [],
          connect2index: [],
        })
      })
      Nmol.enumBonds = {}
      Object.keys(mol.bonds ?? {}).forEach((colorIndex) => {
        // const step = (+100 / mol.bonds[colorIndex].length).toFixed(2)
        Nmol.enumBonds[colorIndex] = reactive({
          array: mol.bonds[colorIndex],
          range: [0, mol.bonds[colorIndex].length],
          keepSame2Index: [],
          connect2index: [],
        })
      })
      console.log('mol:', mol)
      console.log('Nmol:', Nmol)
      return Nmol
    })
  },
)
watch(coresData, () => {
  console.log('coreData', coresData.value)
}, { deep: true })
</script>

<template>
  <div>
    <div v-for="core in coresData" :key="core" class="flex flex-nowarp justify-between m-1">
      <div class="w-25% flex-grow-1 m-1 ml-0">
        <card-rdkit min-w-200px v-bind="core" />
      </div>
      <div class="w-25% justify-around flex-grow-1 place-items-center min-w-420px b-(2 rd-2 solid indigo-2) bg-slate-2 m-1 mr-0">
        <combo-setting :core="core" :ligands="props.ligands" />
      </div>
      <div class="w-70%" />
    </div>
  </div>
</template>

<style scoped>
</style>
