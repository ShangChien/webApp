<script setup lang='ts'>
import { ref, computed, onMounted,watch,nextTick,watchEffect,toRaw,reactive } from "vue"
import { NTree,NSwitch,NSpace,NSlider } from 'naive-ui'
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import classTree from "@/components/enumMole/classTree.vue"
import { useEnumStore } from '@/stores/enumStore'
import type { molData } from "@/components/types"
//use&dataState
const options:any= Array.from(Array(10).keys()).map((i:number)=> {return {label:i,value:i}}) 
const props=defineProps<{ligands:molData[]|any,cores:molData[]|any}>()
const enumStore = useEnumStore()
const coresData=ref()
const ligandData=ref()
const range=ref([100,100])
//控制显示core正在设置的高亮原子
const iAtomColor=ref(0)
//同上。。
const iBondColor=ref(0)

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
	<div v-for="core in coresData" class="flex justify-between m-1" >
		<div class="w-25% min-w-375px flex-grow-1 m-1 ml-0"><card-rdkit v-bind="core" ></card-rdkit></div>
		<div class="w-50% justify-around flex-grow-1 place-items-center relative b-2 rd-2 min-w-480px b-indigo-100 m-1 mr-0">
			<div v-for="(v,k,i) in core.enumAtoms" >
				<div v-if="iAtomColor==k" class="bg-slate-2 b-2 rd-2 b-indigo-50 m-1 min-w-465px ">
					<div class="flex flex-nowrap justify-begin m-1 mb-0 " >
						<div class="font-900 w-120px text-center m-1">atoms组合设置 :</div>
						<div class="flex flex-nowarp offset mt--1">
							<div v-for="option in Object.keys(core.enumAtoms)"
								:style="{color:'hsla('+(+option+8.6)*36+',90%,60%,1)'}" 
								@click="iAtomColor=+option">
								<div v-if="+option == iAtomColor" class="b-2 rd-t-1.5 rd-b-0 b-indigo-2 b-b-slate-50 bg-slate-50">
									<div class="i-carbon-circle-filled text-3xl opacity-100 m-2px"></div>
								</div>
								<div v-else >
									<div class="i-ic-sharp-circle text-3xl opacity-20 m-4px " ></div>
								</div>
							</div>
						</div>
					</div>
					<div class="gridView place-items-center v-full min-w-450px b-2 rd-1.5 b-indigo-2 bg-slate-50 m-1 mt-0">
						<span class="font-400">对接位点:</span>
						<div class="b-2 rd-2 b-indigo-50 bg-white flex flex-nowrap">
							<div v-for="option in options" 
       					:key="option.value" 
       					@click="v.connect2index.includes(option.value) ? 
									v.connect2index.splice(v.connect2index.indexOf(option.value),1) :
									v.connect2index.push(option.value)"
       					:style="{color:'hsla('+(option.value+8.6)*36 + ',90%,60%,1)'}">
								<div v-if="v.connect2index.includes(option.value)" class="i-carbon-circle-filled text-3xl m-1"></div>
								<div v-else class="i-ic-sharp-circle text-3xl opacity-20 m-1
									hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" ></div>
							</div>
						</div>
						<span class="font-400">相同设置:</span>
						<div  class="b-2 rd-2 b-indigo-50 bg-white flex flex-nowrap">
							<div v-for="option in Object.keys(core.enumAtoms)"
								:key="+option" 
								:style="{color:'hsla('+(+option+8.6)*36+',90%,60%,1)'}"
								@click="v.keepSame2Index.includes(option) ? 
									v.keepSame2Index.splice(v.keepSame2Index.indexOf(option),1) :
									v.keepSame2Index.push(option)" >
								<div v-if="v.keepSame2Index.includes(option)" class="i-carbon-circle-filled text-3xl m-1"></div>
								<div v-else class="i-ic-sharp-circle text-3xl opacity-20 m-1
									hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" ></div>
							</div>
						</div>
						<span class="font-400">取代比例:</span>
						<n-slider class="w-90% p-1 b-2 rd-2 b-indigo-50 bg-white"
							v-model:value="v.range" range :format-tooltip="(value: number) => `${value}%`"
							:step="+(100/v.list.length).toFixed(2)" />
					</div>
				</div>
			</div>
			<div v-for="(v,k,i) in core.enumBonds" >
				<div v-if="iBondColor==k" class="bg-slate-2 b-2 rd-2 b-indigo-50 m-1 ">
					<div class="flex flex-nowrap justify-begin m-1 mb-0" >
						<div class="font-900 w-120px text-center m-1">bonds组合设置 :</div>
						<div class="flex flex-nowarp offset mt--1">
							<div v-for="option in Object.keys(core.enumBonds)"
								:style="{color:'hsla('+(+option+8.6)*36+',90%,60%,1)'}" 
								@click="iBondColor=+option">
								<div v-if="+option == iBondColor" class="b-2 rd-t-1.5 rd-b-0 b-indigo-2 b-b-slate-50 bg-slate-50">
									<div class="i-carbon-circle-filled text-3xl opacity-100 m-2px"></div>
								</div>
								<div v-else >
									<div class="i-ic-sharp-circle text-3xl opacity-20 m-4px " ></div>
								</div>
							</div>
						</div>
					</div>
					<div class="gridView place-items-center v-full min-w-450px b-2 rd-1.5 b-indigo-2 bg-slate-50 m-1 mt-0">
						<span class="font-400">对接位点:</span>
						<div class="b-2 rd-2 b-indigo-50 bg-white flex flex-nowrap">
							<div v-for="option in options" 
       					:key="option.value" 
       					@click="v.connect2index.includes(option.value) ? 
									v.connect2index.splice(v.connect2index.indexOf(option.value),1) :
									v.connect2index.push(option.value)"
       					:style="{color:'hsla('+(option.value+8.6)*36 + ',90%,60%,1)'}">
								<div v-if="v.connect2index.includes(option.value)" class="i-carbon-circle-filled text-3xl m-1"></div>
								<div v-else class="i-ic-sharp-circle text-3xl opacity-20 m-1
									hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" ></div>
							</div>
						</div>
						<span class="font-400">相同设置:</span>
						<div  class="b-2 rd-2 b-indigo-50 bg-white flex flex-nowrap">
							<div v-for="option in Object.keys(core.enumBonds)"
								:key="+option" 
								:style="{color:'hsla('+(+option+8.6)*36+',90%,60%,1)'}"
								@click="v.keepSame2Index.includes(option) ? 
									v.keepSame2Index.splice(v.keepSame2Index.indexOf(option),1) :
									v.keepSame2Index.push(option)" >
								<div v-if="v.keepSame2Index.includes(option)" class="i-carbon-circle-filled text-3xl m-1"></div>
								<div v-else class="i-ic-sharp-circle text-3xl opacity-20 m-1
									hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" ></div>
							</div>
						</div>
						<span class="font-400">取代比例:</span>
						<n-slider class="w-90% p-1 b-2 rd-2 b-indigo-50 bg-white"
							v-model:value="v.range" range :format-tooltip="(value: number) => `${value}%`"
							:step="+(100/v.list.length).toFixed(2)" />
					</div>
				</div>
			</div>	
		</div>
		<div class="w-30%"></div>
	</div>
</template>
<style scoped>
.gridView {
  display: grid;
  grid-template-columns: minmax(60px,2fr) minmax(385px,8fr);
  grid-template-rows: 3rem 3rem 3rem ;
}
.gridView1 {
  display: grid;
  grid-template-columns: 80px auto;
  grid-auto-rows: auto;
}
.offset{
	transform: translate(0,2px)
}

</style>