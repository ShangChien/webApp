<script setup lang="ts">
import { refDebounced } from '@vueuse/core'
import { computed, ref, onMounted,onUpdated,onBeforeMount } from "vue";
import type { molData } from "@/components/types";
import cardRdkit from "@/components/rdkitComponent/cardRdkit.vue";
import { NInputNumber,NPagination,NPopover,NTag,NScrollbar } from "naive-ui";
import { useElementSize } from '@vueuse/core'
const props=defineProps<{molList:molData[];cols:number;rows:number}>()
const outBox = ref<HTMLElement | null>(null)
const { width:outBoxW } = useElementSize(outBox)
const initArray = computed(()=>props.molList)
const inputCols = ref(props.cols)
const inputRows = ref(props.rows)
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
<div ref="outBox" class="b-2 rd-2  b-indigo-100 relative min-w-230px min-h-160px">
  <div class="absolute rd-2 z-2 top-0 menubg" :style="{'width':outBoxW+'px'}">
    <div class="ma-1 title-grid">
      <div><n-pagination v-model:page="currentPageInput"
                      :page-count="pageCount"
                      :page-slot="7"
                      :simple='outBoxW < 500 ? true : false'
                      ref="title"
                      class="ma-1 min-w-180px"
                      size="medium"
                      show-quick-jumper>
          <template #goto>跳至:</template>
        </n-pagination></div>
      <div class="justify-self-end mt-1">
        <n-tag v-show='outBoxW > 700 ? true : false' :bordered="false" class=" tagbg  ">
          <span class="text-1.2em ">Total: {{initArray.length}}</span>
          <span class="text-1.2em "> ( {{pageSize}} / page)</span>
        </n-tag>
        <n-popover placement="bottom-end" trigger="click">
          <template #trigger>
            <button class="i-fluent-table-settings-20-filled 
                        text-2.5em
                        c-indigo-300 
                        inline-block
                        mt--1.3
                        hover:c-teal-300
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
  </div>
  <n-scrollbar class="max-h-100vh" >
    <div class="wrapper1 mr-2.5 mt-12 h-100\%" >
      <card-rdkit class="w-100\% h-100\%"
        v-for="(itemInner,indexInner) of arrayShow"
        v-bind="itemInner" 
        :key="indexInner"
      />
    </div>
  </n-scrollbar> 
</div>
</template>
<style scoped>
.wrapper1 {
  display: grid;
  grid-template-columns: repeat(v-bind(cols),1fr);
  grid-template-rows: repeat(v-bind(rows),1fr);
  grid-row-gap: 0.5em;
  grid-column-gap: 0.5em;
  padding:4px;
  padding-top: 0;
}
.menubg {
  border-radius: 6px 6px 0px 0px; 
  backdrop-filter: saturate(70%) blur(12px);
  background: rgb(227,238,255,20%);
}
.tagbg {
  border-radius: 6px 6px 6px 6px; 
  backdrop-filter: saturate(70%) blur(12px);
  background: rgb(227,238,255,10%);
}
/* .box {
  overflow:auto;
  resize:both;
} */
.title-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
}
</style>