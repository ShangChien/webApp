<script setup lang="ts">
import { NTag,NSelect, type SelectRenderTag } from "naive-ui";
import { ref,computed,onMounted,h } from "vue";
import { useVModel } from '@vueuse/core'
import { useEnumStore } from '@/stores/enumStore'
const props=defineProps<{tags:string[]}>()
const emit=defineEmits(['update:tags'])
const tagsValue = useVModel(props,'tags',emit)

const enumStore = useEnumStore()
const tagsOption=computed<{label:string,value:string}[]>(()=>{
	return enumStore.getAllLabels.map((v:string)=>({label:v,value:v}))
})
const renderTag:SelectRenderTag=({ option, handleClose })=>{
	return h(NTag,
        	{
						type:'success',
        	  closable:true,
						round:true,
						size:"small",
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
	
})
</script>
<template>
	<n-select type="info" 
						v-model:value="tagsValue"
						size="small" 
						:options="tagsOption"
						:render-tag="renderTag"
						multiple clearable filterable tag >
	</n-select>
</template>