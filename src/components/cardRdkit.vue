<script setup lang="ts">	
import { onMounted,ref,h,reactive,toRefs,watch } from 'vue';	
import type { Component } from 'vue'
import type { molData } from '@/components/types';
import svgRdkit from '@/components/svgRdkit.vue';
import editableRdkit from '@/components/editableRdkit.vue'
import { NCard,NCheckbox,NSpace,NDropdown,NButton,NIcon,NPopover,NEllipsis,NModal } from 'naive-ui';
import { Dots } from '@vicons/tabler';
import { Edit,Delete } from '@vicons/carbon';
const props = defineProps<molData>()
const emit = defineEmits(['itemChecked','changeItem','itemDeleted'])
//handle Modal
const showModal=ref(false)
function onNegativeClick () {
  showModal.value = false
}
function onPositiveClick () {
  showModal.value = false
}
function handleDropOption(key:string,value:any){
  if (key=='edit'){
    showModal.value = true
  }else if (key=='delete'){
    //
  }else{
    //
  }
}

const show=ref(false)

const renderIcon = (icon: Component) => {
  return () => {
    return h(NIcon, null, {
      default: () => h(icon)
    })
  }
}
const options=[
        {label: 'edit',
          key: 'edit',
          icon: renderIcon(Edit)
        },
        {label: 'delete',
          key: 'delete',
          icon: renderIcon(Delete)
        }
      ]
const checked = ref(false)      
function visible(){
  if (checked.value==true){
    show.value=true
  }
  else{
    show.value=!show.value
  }
  
}
</script>
<template>
<n-card hoverable class="card"
				@mouseover="visible"
				@mouseout="visible">
	<template #cover  >
		<n-space justify="space-between" style="padding-left: 2%;padding-top: 1%;padding-right: 1%;" >
			<n-checkbox v-model:checked="checked" v-show='show' />
      &nbsp
			<div v-show='show'>
        <n-dropdown :options="options" placement="bottom-end" @select="handleDropOption">
          <n-button tertiary circle size="tiny" >
            <template #icon>
              <n-icon><Dots/></n-icon>
            </template>
          </n-button>
			  </n-dropdown>
        <n-modal v-model:show="showModal" :mask-closable="false"> 
          <n-card
              style="width: 800px; height: 600px;"
              title="位点重新选取"
              :bordered="false"
              size="huge"
              role="dialog"
              aria-modal="true">
              <editable-rdkit v-bind="props" qsmiles='*~*' :width="400" :height="300"/>
              <template #footer>
                <n-space>
                  <n-button @Click="onNegativeClick">取消</n-button>
                  <n-button @Click="onPositiveClick">确定</n-button>
                </n-space>
              </template>
          </n-card>
        </n-modal>

      </div> 
	  	</n-space>
      <n-popover trigger="hover" width="trigger" placement="right-end" >
        <template #trigger>
          <svg-rdkit v-bind="props" style="width:100%; height:100%;"/>
        </template>
        <span>
          <n-ellipsis style="max-width: 100%">
            {{props.smiles}}
          </n-ellipsis>
		      atoms:{{props.atoms}}<br>
		      bonds:{{props.bonds}}
        </span>
      </n-popover>
  </template>
</n-card>

</template>
<style>

</style>