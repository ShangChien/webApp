<script setup lang='ts'>
import { ref, computed, onMounted,watch,nextTick,watchEffect,toRaw,reactive } from "vue"
import { NTree,NSwitch,NSpace,NSlider } from 'naive-ui'
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import classTree from "@/components/enumMole/classTree.vue"
import { useEnumStore } from '@/stores/enumStore'
import type { molData } from "@/components/types"
//use&dataState
const Color=(n:number)=>'hsla('+(n+8.6)*36 +',90%,70%,1)'
const options:any= Array.from(Array(10).keys()).map((i:number)=> {return {label:i,value:i}}) 
const props=defineProps<{ligands:molData[]|any,cores:molData[]|any}>()
const enumStore= useEnumStore()
const prop=computed(()=>enumStore.getById(10))
const coresData=ref()
const ligandData=ref()
const range=ref([100,100])

function generateMarks(range:number[]) {
	let ob={100:'100%'}
	console.log('ssds',range[1]/range[0])
	Array.from(Array(Math.round(range[1]/range[0])).keys()).map((_e,i)=>{
		ob[range[0]*i]=range[0]*i+'%'
	})
	return ob
}
watch(
	()=>props.cores,
	(nCores)=>{
		coresData.value = nCores.map((mol:molData)=>{
			const Nmol = JSON.parse(JSON.stringify(mol))
			//必须新建空对象'enumAtoms'，才能深层赋值
			Nmol['enumAtoms']={}
			Object.keys(mol.atoms ?? {}).map((colorIndex)=>{
				const step = (100/mol.atoms[colorIndex].length).toFixed(2)
				console.log(typeof(step))
				const marks = generateMarks([parseInt(step),100])
				Nmol.enumAtoms[colorIndex] = reactive({
					list : mol.atoms[colorIndex],
					range : [step, 100],
					keepSame2Index : [],
					connect2index  : [],
					marks
				})
			})
			Nmol['enumBonds']={}
			Object.keys(mol.bonds ?? {}).map((colorIndex)=>{
				const step=(+100/mol.bonds[colorIndex].length).toFixed(2)
				Nmol.enumBonds[colorIndex] = reactive({
					list : mol.bonds[colorIndex],
					range : [step, 100],
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
watchEffect(()=>{
	console.log(prop.value)
})
</script>
<template>
	<div v-for="core in coresData" class="flex justify-between m-1" >
		<card-rdkit v-bind="core" class="w-20% flex-grow-0"></card-rdkit>
		<div class="flex justify-around flex-grow-1">
			<div v-for="(v,k) in core.enumAtoms" class="flex-grow-1 m-1" >
				<div>
					<div class="i-ic-sharp-circle text-3xl" 
						:style="{color:'hsla('+(+k+8.6)*36+',90%,60%,1)'}"></div>
					<span>对接颜色:</span>
					<div v-for="option in options" 
       			:key="option.value" 
       			@click="v.connect2index.includes(option.value) ? 
							v.connect2index.splice(v.connect2index.indexOf(option.value),1) :
							v.connect2index.push(option.value)"
       			:style="{color:'hsla('+(option.value+8.6)*36 + ',90%,60%,1)'}"
						:class='["site",{active : v.connect2index.includes(option.value)}]'
       			class="i-ic-sharp-circle  text-3xl opacity-20
       			       hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" />
				</div>
				<div>
					<span >取代比例：</span>
					<n-slider class="inline-block" 
						v-model:value="v.range" range 
						:marks="v.marks"
						:step="+(100/v.list.length).toFixed(2)"  />
				</div>
			</div>
		</div>
	</div>
</template>
<style scoped>
.outcircle {
  color: v-bind('Color');
}
.site.active{
	opacity: 1;
}
</style>