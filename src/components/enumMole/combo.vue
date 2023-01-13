<script setup lang='ts'>
import { ref,watch,reactive } from "vue"
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import comboSetting from "@/components/enumMole/comboSetting.vue"
import { useEnumStore } from '@/stores/enumStore'
import type { molData } from "@/components/types"
//use&dataState
const props=defineProps<{ligands:molData[]|any,cores:molData[]|any}>()
const coresData=ref()
watch(
	()=>props.cores,
	(nCores)=>{
		coresData.value = nCores.map((mol:molData)=>{
			const Nmol = JSON.parse(JSON.stringify(mol))
			//必须新建空对象'enumAtoms'，才能深层赋值
			Nmol['enumAtoms']={}
			Object.keys(mol.atoms ?? {}).map((colorIndex)=>{
				const step = (100/mol.atoms[colorIndex].length).toFixed(2)
				Nmol.enumAtoms[colorIndex] = reactive({
					list : mol.atoms[colorIndex],
					range : [+step, 100],
					keepSame2Index : [],
					connect2index  : [],
				})
			})
			Nmol['enumBonds']={}
			Object.keys(mol.bonds ?? {}).map((colorIndex)=>{
				const step=(+100/mol.bonds[colorIndex].length).toFixed(2)
				Nmol.enumBonds[colorIndex] = reactive({
					list : mol.bonds[colorIndex],
					range : [+step, 100],
					keepSame2Index : [],
					connect2index  : [],
				})
			})
			console.log(mol)
			console.log(Nmol)
			return Nmol
		})
	}
)
watch(coresData,()=>{
	console.log(coresData.value)
},{deep:true})

</script>
<template>
<div>
	<div v-for="core in coresData" class="flex flex-nowarp justify-between m-1" >
		<div class="w-25% flex-grow-1 m-1 ml-0"><card-rdkit min-w-200px v-bind="core" ></card-rdkit></div>
		<div class="w-25% justify-around flex-grow-1 place-items-center b-2 rd-2 min-w-420px bg-slate-2 b-indigo-2 m-1 mr-0">
			<combo-setting :dataList="core"></combo-setting>
		</div>
		<div class="w-70%"></div>
	</div>
</div>
</template>
<style scoped>

</style>