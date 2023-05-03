<script setup lang='ts'>
import { ref, computed, onMounted } from "vue"
import axios from "axios";
import type { molData } from "@/components/types";
import { NScrollbar,NPopover,NInputNumber,NSwitch } from "naive-ui";
const props=defineProps<{molList:molData[]}>()
const cdxml4export=ref<string[]>([])
const downloaded=ref<number[]>([])
const cols=ref<number>(10)
const rows=ref<number>(20)
const scale=ref<number>(1)
const singlePage=ref<boolean>(false) 
function dloadCDX(data:string,index:number){
	const now = new Date();
	const fileName = `${index}-${now.getFullYear()}${now.getMonth()+1}${now.getDate()}${now.getHours()}${now.getMinutes()}${now.getSeconds()}.cdxml`;
	const blob = new Blob([data], {type : 'text/plain'});
	const link = document.createElement('a');
	link.href = window.URL.createObjectURL(blob);
	link.download = fileName;
	link.click();
	downloaded.value.push(index)
	window.URL.revokeObjectURL(link.href);
}

function getCdxFile(){
	console.log('POSTing to get cdxfile...',)
	let data =computed(()=>{
		return {
			smiles : props.molList.map((el)=>el.smiles),
			layout : {
				cols  		 :	cols.value,
				rows			 :	rows.value,
				scale 		 :	scale.value,
				singlePage :  singlePage.value
			}
		}
	}) ;
  axios.post(
		'api/cdxml', 
		data.value,
	).then((res:any)=>{
    console.log('receve res')
		cdxml4export.value=res.data.data
  }).catch((error)=>{
    console.log(error);
  });
}

</script>
<template>
	<div class="font-LX-B m--2 mb--1 p--2 ">
    <div class="text-xl">Format:</div>
    <div class="flex flex-nowrap justify-start gap-1.5 m-0 pr-4 p-0">
      <div class="bg-slate-1 rd-2 text-center 
									hover:(bg-slate-2 cursor-pointer)">
        <div class="i-fe:smile text-2xl c-yellow-4" ></div>
        <div class="pl-1 pr-1 ">smiles</div> 
      </div>
      <div class="relative">
				<div class="bg-slate-1 rd-2 text-center 
										hover:(bg-slate-2 cursor-pointer)"
										@click="getCdxFile">
        	<div class="i-mdi:xml text-2xl c-green-4 mb-0 pb-0"  ></div>
        	<div class="pl-1 pr-1 font-LX-B text-neutral-7">cdxml</div>
				</div>
				<div class="absolute rd-5 top--3 right--3 bg-rose-2 
										hover:(bg-rose-3 cursor-pointer)">
					<n-popover placement="bottom-end" trigger="click">
						<template #trigger>
      				<div class='i-ic:outline-settings-suggest text-2xl rd-5 c-rose-6'/>
						</template>
						<template #default>
							<div class="font-LX-B m--2 mb--1 mt--1">
								<div class="text-1.2em bg-slate-1 text-center rd-1 mb-1"
									>cdxml文件布局设置：</div>
								<div class="grid justify-center grid-cols-2 gap-1 w-170px">
              		<span class="self-center justify-self-end">cols:</span>
              		<n-input-number v-model:value="cols" 
              		                :min="1"
              		                size="tiny"
              		                button-placement="both" />
              		<span class="self-center justify-self-end">rows:</span>
              		<n-input-number v-model:value="rows" 
              		                :min="1"
              		                size="tiny"
              		                button-placement="both" />
									<span class="self-center justify-self-end">scale:</span>
              		<n-input-number v-model:value="scale"
																	:min="0" 
              		                size="tiny"
              		                button-placement="both" />
									<span class="self-center justify-self-end">singePage:</span>
              		<n-switch v-model:value="singlePage"
              		                size="small" />			
            		</div>
							</div>
						</template>
					</n-popover>
				</div>
      </div>
    </div>
		<div class="flex flex-nowrap flex-col m-0 p-0 mt-1 rd-2 bg-slate-50">
			<div v-if="cdxml4export.length === 0"
				class="flex flex-col justify-center items-center">
				<div class="i-fxemoji-expressionless text-5xl m-5"></div>
				no data
			</div>
			<div v-else>
				<n-scrollbar class="max-h-40vh p-0 m-0 max-w-full text-xl">
					<div v-for="(item,index) in cdxml4export">
						<div class="flex flex-nowarp justify-around bg-slate-1 m-1 rd-1 hover:cursor-pointer"
							@click="dloadCDX(item,index)">
							<div class="">file-{{index}}</div>
							<div class="i-grommet-icons:download-option m-1"
								:class="downloaded.includes(index) ? ['c-purple-6']:['c-blue-6'] "
							></div>
						</div>
					</div>
				</n-scrollbar>
			</div>
		</div>
  </div>
</template>
<style>
</style>