<script setup lang='ts'>
import { computed, reactive, ref } from 'vue'
import { NButton, NDataTable } from 'naive-ui'
import type { dataResults } from '@/components/types'

const props = defineProps<{ data: dataResults[] }>()

const dom = ref()
const columns: any = computed(() => {
  const cols = Object.keys(props.data[0]).filter(key => key !== 'key')
  console.log(cols)
  return cols.map((e) => {
    if (e === 'name') {
      return { title: 'Name', key: 'name', defaultSortOrder: 'ascend', sorter: 'default', width: 200 }
    } else {
      return { title: e, key: e, sorter: (row1, row2) => row1[e] - row2[e], resizable: true, minWidth: 200 }
    }
  })
})
// const columns: any = [
//   { title: 'Name', key: 'name', defaultSortOrder: 'ascend', sorter: 'default', width: 100 },
//   { title: 'homo', key: 'homo', sorter: (row1, row2) => row1.age - row2.age, resizable: true, minWidth: 100 },
//   { title: 'lumo', key: 'lumo', sorter: (row1, row2) => row1.age - row2.age, resizable: true, minWidth: 100 },
//   { title: 'eg', key: 'eg', sorter: (row1, row2) => row1.age - row2.age, resizable: true, minWidth: 100 },
// ]

const paginationReactive = reactive({
  page: 2,
  pageSize: 9,
  showSizePicker: true,
  pageSizes: [5, 9, 13, 17, 21, 25, 50, 100],
  onChange: (page: number) => {
    paginationReactive.page = page
  },
  onUpdatePageSize: (pageSize: number) => {
    paginationReactive.pageSize = pageSize
    paginationReactive.page = 1
  },
})

function exportCSV(data: dataResults[]) {
  const now = new Date()
  const fileName = `oled-ml-${now.getFullYear()}${now.getMonth() + 1}${now.getDate()}${now.getHours()}${now.getMinutes()}${now.getSeconds()}.csv`
  const columns = Object.keys(data[0])

  // 将列名添加到CSV字符串
  const csv = [
    columns.join(','), // 添加列名行
    ...data.map(obj => columns.map(col => Array.isArray(obj[col]) ? obj[col].join(',') : obj[col]).join(',')), // 添加数据行
  ].join('\n')

  // 创建并下载CSV文件
  const csvData = new Blob([csv], { type: 'text/csv' })
  const csvUrl = URL.createObjectURL(csvData)

  const link = document.createElement('a')
  link.href = csvUrl
  link.download = fileName
  document.body.appendChild(link)
  link.click()
  URL.revokeObjectURL(csvUrl)
  document.body.removeChild(link)
}
</script>

<template>
  <div class="h-full w-full box-border relative flex">
    <NDataTable
      v-if="props.data.length > 0"
      ref="dom"
      class="h-full"
      flex-height
      size="medium"
      :columns="columns"
      :data="props.data"
      :pagination="paginationReactive"
    />
    <div v-else class="ma text-slate text-center">
      <span>mechine learning results</span>
      <div class="i-carbon-face-dizzy text-4xl ma" />
      <span>no data</span>
    </div>
    <NButton :disabled="!(props.data?.length > 0)" class="absolute bottom-1 left-1 " tertiary type="primary" size="small" @click="exportCSV(props.data)">
      export results
    </NButton>
  </div>
</template>

<style>
</style>
