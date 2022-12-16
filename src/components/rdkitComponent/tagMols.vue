<script setup lang="ts">
import { NTag,NSelect, type SelectRenderTag } from "naive-ui";
import { inject,ref,computed,onMounted,h } from "vue";
import { useEnumStore } from '@/stores/enumStore'
const enumStore = useEnumStore()
const tagSelected =ref ([])
const tagsOption=computed<{label:string,value:string}[]>(()=>enumStore.getAllLabels.map((v:string)=>({label:v,value:v})))
//const labels =computed<{label:string,value:string}[]>(()=>enumStore.getAllLabels.map((v:string)=>({label:v,value:v})))
//const tags=ref([])
const labels=enumStore.getAllLabels
//tagsOption.value=labels.map((v:string,i:number)=>({label:i+1+' '+v,value:v}))
const renderTag:SelectRenderTag=({ option, handleClose })=>{
	return h(NTag,
        	{
						type:'success',
        	  closable:true,
        	  onMousedown: (e:FocusEvent)=>{ 
							e.preventDefault() 
						},
					 	onClose: (e:MouseEvent)=>{
        	    e.stopPropagation()
        	    handleClose()
        	  }
					},
        	{ default: () => option.label } )
}
onMounted(()=>{
	
	console.log(labels,tagsOption.value)
})
</script>
<template>
	<n-select type="info" 
						v-model="tagSelected"
						size="small" 
						:options="tagsOption"
						:render-tag="renderTag"
						multiple clearable filterable tag >
	</n-select>
</template>