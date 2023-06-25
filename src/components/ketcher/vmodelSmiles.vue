<script setup lang='ts'>
import { ref, computed, onMounted } from "vue"
import { NButton,NModal,NCard } from "naive-ui"
import initKetcher from './initKetcher.vue'

const smiles = defineModel<string>('smiles')
const showModal = ref(false)

</script>
<template>
<div class="flex-none flex flex-nowrap justify-between items-center
h-24px m-0 p-0 rd-1 b-blue-2 box-border b-2
active:(outline outline-2px outline-blue-2)
bg-slate-2 text-gray-8
hover:bg-blue-1" >
  <input ref="input" class='flex-auto text-lg outline-0 b-0 rd-1 pl-1 m-0 p-0 box-border h-full'
    v-model="smiles" placeholder="Add conditions by button below"/>
  <div v-show="!!smiles" class="aspect-ratio-1 rd-50% h-5 mr-1 bg-slate-100
  flex justify-center items-center hover:(cursor-pointer bg-gray-200)"
  @click="(e:any)=>{smiles='';e.target.previousElementSibling.focus()}">
    <div class="i-ion-close"></div>
  </div>
  <div @click="showModal=true"
    class="aspect-ratio-1 rd-50% h-7 box-border m-1 mr-2 bg-slate-100 
      flex justify-center items-center hover:(cursor-pointer bg-gray-200)">
    <div class="i-carbon-cloud-satellite text-xl bg-rose-500"></div>
  </div>
  <n-modal v-model:show="showModal" display-directive="show">
    <n-card class="w-80vw h-80vh relative" 
      :bordered = "false"
      size="huge"
      role="dialog"
      aria-modal="true">
      <template #default>
        <init-ketcher class="w-full h-full"
          v-model:smiles="smiles"
        ></init-ketcher>
      </template>
      <template #footer>
        <div class="flex justify-center gap-2">
          <n-button @Click="showModal = false">取消</n-button>
          <n-button @Click="showModal = false">确定</n-button>
        </div>
      </template>
    </n-card>
  </n-modal>
</div>
</template>
<style>
</style>