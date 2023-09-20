<script setup lang='ts'>
import { ref } from 'vue'
import nglViewer from '@/components/nglViewer.vue'
import mlFlow from '@/components/ml/mlFlow.vue'

const molStr = ref<string>('')
const result = ref<{ name: string; homo: number; lumo: number; eg: number }[]>([])
</script>

<template>
  <div class="flex flex-col justify-start items-center w-full">
    <div class="flex-none flex flex-nowrap justify-between items-center w-full">
      <nglViewer :data="molStr" />
      <div>
        <div>
          <div>predict all</div>
          <div>predict current</div>
        </div>
        <div>
          <n-checkbox-group v-model:value="cities">
            <n-space item-style="display: flex;">
              <n-checkbox value="Beijing" label="北京" />
              <n-checkbox value="Shanghai" label="上海" />
              <n-checkbox value="Guangzhou" label="广州" />
              <n-checkbox value="Shenzen" label="深圳" />
            </n-space>
          </n-checkbox-group>
        </div>
      </div>
    </div>
    <mlFlow class="flex-auto" :data="result" />
  </div>
</template>

<style>
</style>
