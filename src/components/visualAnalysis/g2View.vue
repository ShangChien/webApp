<script setup lang='ts'>
import { onMounted, ref } from 'vue'
import { Chart } from '@antv/g2'

const g2View = ref()
onMounted(() => {
  const chart = new Chart({
    container: g2View.value,
    autoFit: true,
  })

  const data = [
    { genre: 'Sports', sold: 275 },
    { genre: 'Strategy', sold: 115 },
    { genre: 'Action', sold: 120 },
    { genre: 'Shooter', sold: 350 },
    { genre: 'Other', sold: 150 },
  ]

  chart
    .interval()
    .data(data)
    .encode('x', 'genre')
    .encode('y', 'sold')
    .encode('key', 'genre')
    .animate('update', { duration: 300 })

  // 渲染可视化
  chart.render()
})
</script>

<template>
  <div ref="g2View" />
</template>

<style>
</style>
