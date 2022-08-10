<script setup lang="ts">
import { NIcon, NButton, NGradientText } from "naive-ui";
import { Theater } from "@vicons/carbon";
import { ref } from "vue";

const siteType = ref(0)
const Color = (n:number) => {
  const colorAngle = Math.floor((n+4.6)*36);
  return 'hsla('+ colorAngle +',90%,70%,1)';
}
</script>
<template>
<div style="margin-top:-5px">
	<n-gradient-text 
		style="font-weight:600; background-image: linear-gradient(80deg, #414f88, #7a6176, #a57663, #cc8b4b)">
    Type:
  </n-gradient-text>
	<n-button circle  quaternary
		:style="{left:5+'px',top:4+'px', background:Color(siteType)}"
		@click="siteType++" >
    <template #icon>
      <n-icon size="20" color="#0e7a0d" ><theater /></n-icon>
    </template>
  </n-button>
</div>
</template>
<style>
</style>