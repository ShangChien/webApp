<script setup lang="ts">
import type { molData } from "@/components/types";
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import { computed, onMounted } from "vue";
const props = defineProps<{molList:molData[]}>();
const cols=computed(()=>props.molList.length)
onMounted(()=>{
  console.log(cols.value)
})
</script>
<template>
<div class="wrapper2"  >
  <card-rdkit 
    v-for="(item,index) of props.molList"
    v-bind="item" 
    :key="index"
  />
</div>
</template>
<style scoped>
.wrapper2 {
  display: grid;
  grid-template-columns: repeat( v-bind('cols'),1fr);
  grid-column-gap: 0.5em;
}
</style>
