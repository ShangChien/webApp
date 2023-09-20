<script setup lang='ts'>
import { ref } from 'vue'
import { NDataTable } from 'naive-ui'

const props = defineProps<{ data: { name: string; homo: number; lumo: number; eg: number }[] }>()

const dom = ref()
const columns = [
  { title: 'Name', key: 'name', defaultSortOrder: 'ascend', sorter: 'default' },
  { title: 'homo', key: 'homo', sorter: (row1, row2) => row1.age - row2.age },
  { title: 'lumo', key: 'lumo', sorter: (row1, row2) => row1.age - row2.age },
  { title: 'eg', key: 'eg', sorter: (row1, row2) => row1.age - row2.age },
]
</script>

<template>
  <div class="h-full w-full">
    <NDataTable
      v-if="props.data.length > 0"
      ref="dom"
      :columns="columns"
      :data="props.data"
      :pagination="{ pageSize: 5 }"
    />
    <div v-else class="flex justify-center items-center">
      <div class="i-noto-hamburger animate-bounce-alt animate-count-infinite animate-duration-1.5s text-4xl" />
      <span>no data</span>
    </div>
  </div>
</template>

<style>
</style>
