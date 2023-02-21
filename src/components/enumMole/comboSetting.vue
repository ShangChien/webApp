<script setup lang='ts'>
import { ref, computed, onMounted,inject } from "vue"
import { NButton } from "naive-ui";
import axios from "axios";
import comboSettingItem from "@/components/enumMole/comboSettingItem.vue"
import type { enumData,molData } from "@/components/types"
//use&dataState
const props=defineProps<{
	core:{
		enumAtoms:enumData,
		enumBonds:enumData,
	},
	ligands:molData[]
}>()
const setting=['enumAtoms','enumBonds']
const settingName=ref('enumAtoms')

//axios请求后端
const enumResult:any = inject('enumResult')
function enumMol(){
	console.log('axios',props)
	//axios.get('api/todos/1').then(json => console.log(json))
  axios.post(
		'api/enum', 
		props,
		// 
	).then(function (response:any) {
		let data=response.data.data
    enumResult.value=data.map((el,index)=>{
			return {
				id:index,
				smiles:el
			}
		});
    console.log(response)
  }).catch(function (error) {
    console.log(error);
  });
}

</script>
<template>
	<div class="flex-col flex h-full">
		<div class="flex-none flex flex-nowrap justify-between m-1 mb-0  text-center" :style="{transform:'translate(0,2px)'}">
			<div class="flex-none flex flex-nowrap justify-begin">
				<div class="font-600 text-1.2em p-1 inline-block ">组合设置 :</div>
				<div class="flex flex-nowrap justify-begin">
					<div v-for="item in setting">
						<div class="p-1 rd-t-2 b-2 b-b-0 text-1.2em"
								 :class="[settingName===item ? 'bg-slate-1 b-sky-2' : 'b-slate-2/0']"
								 @click="settingName=item">{{item.substring(4)}}
						</div>
					</div>
				</div>
			</div>
			<n-button type="info" size="small" @click="enumMol">枚举</n-button>
		</div>
		<combo-setting-item v-for="item in setting" 
			class="flex-auto mb-2" 
			v-show="settingName===item" 
			:data="props.core[item]"></combo-setting-item>
	</div>
</template>
<style>
</style>