<script setup lang="ts">	
import { onMounted,ref,h } from 'vue';	
import type { Component } from 'vue'
import type { molData } from '@/components/types';
import svgRdkit from '@/components/svgRdkit.vue';
import { NCard,NCheckbox,NSpace,NDropdown,NButton,NIcon } from 'naive-ui';
import { Dots } from '@vicons/tabler';
import { Edit,Delete } from '@vicons/carbon';
const props = defineProps<molData>()
const show=ref(false)
const checked = ref(false)
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
onMounted(()=>{
	show.value=false
})
</script>
<template>
<n-card hoverable
				@mouseenter="show=true"
				@mouseleave="show=false">
	<template #cover class="svg"  >
		<n-space justify="space-between" style="padding-inline: 2%" >
			<n-checkbox v-model="checked" />
			<n-dropdown :options="options" >
				<n-button text style="font-size: 20px">
					<n-icon><dots/></n-icon>
				</n-button>
			</n-dropdown>
		</n-space>
    <svg-rdkit v-bind="props" style="width:100%; height:100%;"/>
		atoms:{{props.atoms}}<br>
		bonds:{{props.bonds}}
  </template>
</n-card>
</template>
<style>

</style>