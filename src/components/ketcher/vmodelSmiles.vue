<script setup lang='ts'>
import { ref } from 'vue'
import { NButton, NCard, NModal } from 'naive-ui'
import initKetcher from './initKetcher.vue'

const smiles = defineModel<string>('smiles')
const showModal = ref(false)
</script>

<template>
  <div
    class="flex-none flex flex-nowrap justify-between items-center
h-24px m-0 p-0 rd-1 b-blue-2 box-border b-2
active:(outline outline-2px outline-blue-2)
bg-slate-2 text-gray-8
hover:bg-blue-1"
  >
    <input
      v-model="smiles"
      class="flex-auto text-lg outline-0 b-0 rd-1 pl-1 m-0 p-0 box-border h-full" placeholder="Add conditions by button below"
    >
    <div
      v-show="!!smiles" class="aspect-ratio-1 rd-50% h-5 mr-1 bg-slate-100
  flex justify-center items-center hover:(cursor-pointer bg-gray-200)"
      @click="(e:any) => { smiles = '';e.target.previousElementSibling.focus() }"
    >
      <div class="i-ion-close" />
    </div>
    <div
      class="aspect-ratio-1 rd-50% h-7 box-border m-1 mr-2 bg-slate-100
      flex justify-center items-center hover:(cursor-pointer bg-gray-200)"
      @click="showModal = true"
    >
      <div class="i-carbon-cloud-satellite text-xl bg-rose-500" />
    </div>
    <NModal v-model:show="showModal" display-directive="show">
      <NCard
        class="w-80vw h-80vh relative"
        :bordered="false"
        size="huge"
        role="dialog"
        aria-modal="true"
      >
        <template #default>
          <init-ketcher
            v-model:smiles="smiles"
            class="w-full h-full"
          />
        </template>
        <template #footer>
          <div class="flex justify-center gap-2">
            <NButton @click="showModal = false">
              取消
            </NButton>
            <NButton @click="showModal = false">
              确定
            </NButton>
          </div>
        </template>
      </NCard>
    </NModal>
  </div>
</template>

<style>
</style>
