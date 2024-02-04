<script setup lang='ts'>
import { onMounted, onUnmounted, ref, watch } from 'vue'
import { useElementSize } from '@vueuse/core'
import { Chart } from '@antv/g2'
import type { ViewType } from './types'

const props = defineProps<{ option: ViewType }>()
const container = ref()
const { width, height } = useElementSize(container)
let chart: Chart

function initChart(container) {
  const chart = new Chart({ container, autoFit: true })
  return chart
}
function updateChart(chart: Chart, config: ViewType): Promise<void> {
  return new Promise((resolve, reject) => {
    try {
      chart.options(config as any)
      chart.render()
      resolve()
    } catch (error) {
      reject(error)
    }
  })
}

onMounted(() => {
  chart = initChart(container.value)
  watch([() => props.option, width, height], () => {
    updateChart(chart, props.option)
      .catch(error => console.error(props.option, error))
  }, {
    immediate: true,
    flush: 'post',
  })
  // console.log(props.option)
})

onUnmounted(() => chart.destroy())
</script>

<template>
  <div ref="container" class="flex-auto h-full w-full" />
</template>

<style>
</style>
