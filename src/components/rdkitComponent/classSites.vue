<script setup lang="ts">
import { NIcon, NButton, NGradientText,NPopselect } from "naive-ui";
import { Theater } from "@vicons/carbon";
import { h,inject } from "vue";

const siteType:any= inject('siteType')
const options:any= Array.from(Array(10).keys()).map((i:number)=> {return {label:i,value:i}})
                
const Color = (n:any) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,70%,1)'
const renderItem=( option:any )=>{
  return [h(
            NButton,
            { 
              size:'small',
              circle:true, 
              quaternary:true,
	  	        style:{background:Color(option.value)}
            },
            {
              icon: () => {
                return h(
                  NIcon,
                  {
                    style:{size:20,color:"#0e7a0d"}
                  },
                  {
                   default:()=>h( Theater)
                  }
                )}
            }
          )]
}
</script>
<template>
<div style="margin-top:-5px">
	<n-gradient-text 
		style="font-weight:600; background-image: linear-gradient(80deg, #414f88, #7a6176, #a57663, #cc8b4b)">
    Type:
  </n-gradient-text>
  <n-popselect
    v-model:value="siteType"
    :overlap="true"
    placement="bottom-start"
    :render-label="renderItem"
    :options="options"
    size="medium"
    trigger="click"
    scrollable
  >  
  <n-button circle quaternary size='small'
	  	:style="{left:5+'px',top:4+'px', background:Color(siteType)}"
	  	>
      <template #icon>
        <n-icon size="20" color="#0e7a0d" ><theater /></n-icon>
      </template>
    </n-button>
  </n-popselect>
</div>
</template>
