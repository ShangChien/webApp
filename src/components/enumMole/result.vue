<script setup lang='ts'>
import { ref, computed, onMounted } from "vue"
import gridPage from "@/components/rdkitComponent/gridPage.vue";
import { NButton,NScrollbar } from 'naive-ui'
import axios from "axios";
import type { molData } from "@/components/types"
import type { Ref } from 'vue'
import { reactive } from 'vue'
import { useFileSystemAccess } from '@vueuse/core'
//use&dataState
const props = defineProps<{mols:molData[]}>()
//grid layout
const cols=ref<number>(10)
const rows=ref<number>(3)
const cdxFile=ref()
const dataType = ref('Text') as Ref<'Text' | 'ArrayBuffer' | 'Blob'>
const res = useFileSystemAccess({
  dataType,
  types: [{
    description: 'text',
    accept: {
      'text/plain': ['.txt', '.html'],
    },
  }],
  excludeAcceptAllOption: true,
})
const content = res.data
const str =JSON.stringify(reactive({
  isSupported: res.isSupported,
  file: res.file,
  fileName: res.fileName,
  fileMIME: res.fileMIME,
  fileSize: res.fileSize,
  fileLastModified: res.fileLastModified,
}))
async function onSave() {
  await res.save()
}
function download(data:string){
	const blob = new Blob([data], {type : 'text/plain'});
	const fileName = `${new Date().valueOf()}.cdxml`;
	const link = document.createElement('a');
	link.href = window.URL.createObjectURL(blob);
	link.download = fileName;
	link.click();
	window.URL.revokeObjectURL(link.href);
}

function getCdxFile(){
	console.log('axios',props.mols)
	let data = props.mols.map((el)=>el.smiles);
  axios.post(
		'api/cdxml', 
		data,
	).then((res:any)=>{
    console.log('get res',res)
		download(data=res.data.data[0])
  }).catch(function (error) {
    console.log(error);
  });
}


</script>
<template>
<!-- 
    <div class="flex-none flex flex-nowrap justify-between">
		  <div class="flex-none flex flex-nowrap justify-start ">
		  	<n-button size="small" class="rd-1 m-1 bg-slate-2" @click="res.open()">
          Open
        </n-button>
        <n-button size="small" class="rd-1 m-1 bg-slate-2" @click="res.create()">
          New file
        </n-button>
        <n-button size="small" class="rd-1 m-1 bg-slate-2" :disabled="!res.file.value" @click="onSave">
          Save
        </n-button>
        <n-button size="small" class="rd-1 m-1 bg-slate-2" :disabled="!res.file.value" @click="res.saveAs()">
          Save as
        </n-button>
		  </div>
      <div>
        <n-button size="small" class="rd-1 m-1 bg-slate-2"  @click="getCdxFile">
          download result
        </n-button>
      </div>
    </div> -->
  <div class="mt-2">
	  <grid-page class="h-72vh"
	  	:molList="props.mols" 
	  	:cols='cols'
	  	:rows='rows'
      :maxH="70" />
  </div>
</template>
<style>
</style>