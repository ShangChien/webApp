<script setup lang='ts'>
import { ref, computed, onMounted } from "vue"
import comboSettingItem from "@/components/enumMole/comboSettingItem.vue"
import type { enumData } from "@/components/types"
//use&dataState
const props=defineProps<{dataList:{
	enumAtoms:{enumD:enumData},
	enumBonds:{enumD:enumData}
}}>()
const setting=['enumAtoms','enumBonds']
const settingName=ref('enumAtoms')
</script>
<template>
	<div class="flex-col flex h-full">
		<div class="flex-none flex flex-nowrap justify-begin m-1 mb-0 ml-2 text-center" :style="{transform:'translate(0,2px)'}">
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
		<combo-setting-item v-for="item in setting" 
			class="flex-auto mb-2" 
			v-show="settingName===item" 
			:enumD="props.dataList[item]"></combo-setting-item>
	</div>
</template>
<style>
</style>