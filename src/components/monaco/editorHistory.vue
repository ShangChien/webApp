<script setup lang='ts'>
import type { Ref } from 'vue'
import dayjs from 'dayjs'

const props = defineProps<{ records: { text: string; stamp: number; detial: string } | any }>()
const emits = defineEmits<{
  backTo: [index:number]
  diffTo: [index:number]
}>()
const dateStr = (stamp: number) => dayjs(stamp).format('YY.MM.DD HH:mm:ss')
</script>

<template>
  <div>
    <div v-for="(i, k) in props.records" :key="k">
      {{ i.timestamp }}
    </div>
  </div>
</template>

<style>
</style>
