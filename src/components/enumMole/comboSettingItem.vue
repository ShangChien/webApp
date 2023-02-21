<script setup lang='ts'>
import { ref,computed,watchEffect } from "vue"
import { NSlider } from 'naive-ui'
import type { enumData } from "@/components/types"

const props=defineProps<{data:enumData}>()
const keys=computed(()=>Object.keys(props.data))
const options:any= Array.from(Array(10).keys()).map((i:number)=> {return {label:i,value:i}}) 
//初始化显示的索引
const indexN=ref()
watchEffect(()=>{
	indexN.value = (+keys.value[0])
})
</script>
<template>
<div>
	<div v-for="(v,k,i) in props.data" :key="i" 
		:class="[indexN==k ? 'h-full':'h-0']" >
		<div v-if="indexN==k" 
			class="bg-slate-1 rd-2 b-2 b-sky-2 m-1 mt-0 flex-col flex"
			:class="[indexN==k ? 'h-full':'h-0']">
			<div class="flex flex-nowarp mt-1 ml-3">
				<div v-for="option in keys"
					:style="{color:'hsla('+(+option+8.6)*36+',90%,60%,1)'}" 
					@click="indexN=+option">
					<div v-if="+option == indexN" class="rd-t-1.5 rd-b-0 bg-white">
						<div class="i-carbon-circle-filled text-3xl opacity-100 m-1"></div>
					</div>
					<div v-else >
						<div class="i-ic-sharp-circle text-3xl opacity-50 m-1"></div>
					</div>
				</div>
			</div>
			<div class="flex-auto gridView place-items-center min-w-400px rd-1.5 bg-white m-1 mt-0">
				<span class="font-400">对接位点:</span>
				<div class="b-2 rd-2 b-indigo-50 bg-slate-50 flex flex-nowrap mr-2">
					<div v-for="option in options" 
		 				:key="option.value" 
		 				@click="v.connect2index.includes(option.value) ? 
							v.connect2index.splice(v.connect2index.indexOf(option.value),1) :
							v.connect2index.push(option.value)"
		 				:style="{color:'hsla('+(option.value+8.6)*36 + ',90%,60%,1)'}">
						<div v-if="v.connect2index.includes(option.value)" class="i-carbon-circle-filled text-3xl m-0.2"></div>
						<div v-else class="i-ic-sharp-circle text-3xl opacity-20 m-0.2
							hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" ></div>
					</div>
				</div>
				<span class="font-400">相同设置:</span>
				<div  class="b-2 rd-2 b-indigo-50 bg-slate-50 flex flex-nowrap">
					<div v-for="option in keys"
						:key="+option" 
						:style="{color:'hsla('+(+option+8.6)*36+',90%,60%,1)'}"
						@click="v.keepSame2Index.includes(+option) ? 
							v.keepSame2Index.splice(v.keepSame2Index.indexOf(+option),1) :
							v.keepSame2Index.push(+option)" >
						<div v-if="v.keepSame2Index.includes(+option)" class="i-carbon-circle-filled text-3xl m-0.2"></div>
						<div v-else class="i-ic-sharp-circle text-3xl opacity-20 m-0.2
							hover:(i-carbon-circle-filled duration-210 opacity-70 ease-in-out)" ></div>
					</div>
				</div>
				<span class="font-400">取代比例:</span>
				<n-slider class="w-90% p-1 b-2 rd-2 b-indigo-50 bg-slate-50"
					v-model:value="v.range" range :max="v.list.length" 
					:format-tooltip="(value: number) => `${(value/v.list.length)*100}%`"
					:step="1" />
			</div>
		</div>
	</div>
	<div v-if="keys.length==0" class="bg-slate-1 rd-2 b-2 b-sky-2 m-1 mt-0 h-full flex flex-col justify-center items-center">
		<div class="i-fxemoji-expressionless text-5xl  "></div>
		<span>没有标记相应位点</span>
	</div>
</div>	
</template>
<style scoped>
.gridView {
  display: grid;
  grid-template-columns: minmax(60px,2fr) minmax(300px,8fr);
  grid-template-rows: 1fr 1fr 1fr ;
}
</style>