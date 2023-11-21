<script setup lang='ts'>
import { computed, onMounted, onUnmounted, ref } from 'vue'
import { Chart } from '@antv/g2'
import viewSetting from './viewSetting.vue'
import type { ViewType } from './types'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
}>()
const fileIndex = defineModel<number | null>('index', { default: 0 })
const fileOptions = computed(() => props.allFiles.map((el, index) => ({ label: el.name, value: index })))
const container = ref()
let chart: Chart

const config = ref<ViewType>()

function initChart(container) {
  const chart = new Chart({ container, autoFit: true })
  return chart
}

function updateChart(chart: Chart, config: ViewType) {
  chart.options(config)
  chart.render()
}

onMounted(() => {
  chart = initChart(container.value)
  updateChart(chart, config.value)
  console.log(fileOptions)
})

onUnmounted(() => chart.destroy())
</script>

<template>
  <div class="h-full w-full box-border flex justify-center items-center">
    <div class="box-border w-70% h-full flex-auto b-2 rd-1 b-indigo-100 b-solid">
      <div ref="container" class="flex h-full w-full" />
    </div>
    <div class="flex-auto bg-gray-50 h-full w-30%  box-border">
      <viewSetting :all-files="allFiles" :file-index="fileIndex" :config="config" />
    </div>
  </div>
</template>

<style>
</style>
