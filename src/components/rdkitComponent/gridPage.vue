<script setup lang="ts">
import { refDebounced } from '@vueuse/core'
import { computed, ref, onMounted,onUpdated,onBeforeMount } from "vue";
import type { molData } from "@/components/types";
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import { NInputNumber,NPagination,NPopover,NTag } from "naive-ui";

const props=defineProps<{mollist:molData[]}>()
const initArray=computed(()=>props.mollist)
const inputCols = ref(8)
const inputRows = ref(8)
const cols = refDebounced(inputCols, 500)
const rows = refDebounced(inputRows, 500)
const pageSize = computed(()=>cols.value*rows.value)
const pageCount = computed(()=>Math.ceil(initArray.value.length/pageSize.value))
const currentPageInput=ref(1)
const currentPage = refDebounced(currentPageInput, 500)
const arrayShow=computed(()=>{
  let start = (currentPage.value-1)*pageSize.value
  let end = currentPage.value*pageSize.value
  return initArray.value.slice(start,end)
})

</script>

<template>
<div class="b-2 rd-2 b-indigo-1 ">
  <div class="grid grid-cols-2 ">
    <div>
      <n-pagination v-model:page="currentPageInput"
                    :page-count="pageCount"
                    :page-slot="7"
                    class="ma-1 "
                    size="medium"
                    show-quick-jumper>
        <template #goto>跳至:</template>
      </n-pagination>
    </div>
    <div class="justify-self-end ma-1">
      <n-tag :bordered="false" class="mr-1 bg-red-50">
        <span class="text-1.2em ">Total: {{initArray.length}}</span>
      </n-tag>
      <n-tag :bordered="false" class="mr-1 bg-red-50">
        <span class="text-1.2em "> ( {{pageSize}} / page)</span>
      </n-tag>
      <n-popover trigger="click">
        <template #trigger>
          <button class="i-fluent-table-settings-20-filled 
                      text-2.5em
                      c-indigo-400 
                      mt--1.2
                      hover:c-teal-400
                      focus:(c-teal-400 transform-scale-110)" />
        </template>
        <template #default>
          <div class="grid justify-items-end grid-cols-2 gap-2">
            <span class=" pr-1 pl-2 text-xl">cols:</span>
            <n-input-number v-model:value="inputCols" 
                            :update-value-on-input="false"
                            :min="2"
                            size="small"
                            class="w-20 " 
                            button-placement="both" />
            <span class=" pr-1 pl-2 text-xl ">rows:</span>
            <n-input-number v-model:value="inputRows" 
                            :update-value-on-input="false"
                            :min="2"
                            size="small"
                            class="w-20" 
                            button-placement="both" />
          </div>
        </template>
      </n-popover>
    </div>
  </div>
  <div class="wrapper1" >
      <card-rdkit class="w-100% h-100%"
        v-for="(itemInner,indexInner) of arrayShow"
        v-bind="itemInner" 
        :key="indexInner"
      />
  </div>
</div>
</template>
<style scoped>
.wrapper1 {
  display: grid;
  grid-template-columns: repeat( v-bind('cols'),1fr);
  grid-template-rows: repeat(v-bind(rows),1fr);
  grid-row-gap: 0.5em;
  grid-column-gap: 0.5em;
  padding:4px;
  padding-top: 0;
}
</style>