<script setup lang='ts'>
import { reactive, ref } from 'vue'
import { NButton, NDataTable } from 'naive-ui'

const props = defineProps<{ data: { name: string; homo: number; lumo: number; eg: number }[] }>()

const dom = ref()
const columns: any = [
  { title: 'Name', key: 'name', defaultSortOrder: 'ascend', sorter: 'default', width: 100 },
  { title: 'homo', key: 'homo', sorter: (row1, row2) => row1.age - row2.age, resizable: true, minWidth: 100 },
  { title: 'lumo', key: 'lumo', sorter: (row1, row2) => row1.age - row2.age, resizable: true, minWidth: 100 },
  { title: 'eg', key: 'eg', sorter: (row1, row2) => row1.age - row2.age, resizable: true, minWidth: 100 },
]

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

function exportCSV(data: { name: string; homo: number; lumo: number; eg: number }[]) {
  const now = new Date()
  const fileName = `oled-ml-${now.getFullYear()}${now.getMonth() + 1}${now.getDate()}${now.getHours()}${now.getMinutes()}${now.getSeconds()}.csv`
  const columns = Object.keys(data[0])

  // 将列名添加到CSV字符串
  const csv = [
    columns.join(','), // 添加列名行
    ...data.map(obj => columns.map(col => obj[col]).join(',')), // 添加数据行
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
  <div class="h-full w-full box-border relative">
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
    <div v-else class="flex justify-center items-center">
      <div class="i-noto-hamburger animate-bounce-alt animate-count-infinite animate-duration-1.5s text-4xl" />
      <span>no data</span>
    </div>
    <NButton :disabled="!(props.data?.length > 0)" class="absolute bottom-1 left-1 " tertiary type="primary" size="small" @click="exportCSV(props.data)">
      export results
    </NButton>
  </div>
</template>

<style>
</style>
