<script setup lang='ts'>
import { computed, defineAsyncComponent, h, ref, watch } from 'vue'
import { NButton, NSelect, NTag } from 'naive-ui'
import type { Spectrum } from './dataCheckStore'
import { file2Data, normalize, tickInfo } from '@/components/visualAnalysis/utils'
import figCheckOption from '@/components/dataCheck/figCheckOption.vue'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
  index: number
}>()

const g2View = defineAsyncComponent(() => import('@/components/visualAnalysis/g2View.vue'))

const dataCSV = computed(() => {
  let data = file2Data<object>(props.allFiles[props.index])
  const cols = Object.keys(data[0])
  data = normalize(data)
  return { data, cols }
})

// xcol, ycol选项
const Xcol = ref(0)
const XcolName = computed(() => dataCSV.value.cols[Xcol.value])
const Ycol = ref(1)
const YcolName = computed(() => dataCSV.value.cols[Ycol.value])
const colsOptions = computed(() => {
  const cols = dataCSV.value.cols.map((v, i) => ({ label: v, value: i, type: 'info' }))
  return cols
})
function renderTag({ option }) {
  return h(
    NTag,
    {
      type: option.type,
      size: 'small',
    },
    { default: () => option.label },
  )
}

const data4ref = ref<Spectrum[]>([])
const data4check = ref<Spectrum[]>([])
const data4show = computed(() => {
  return [...data4check.value, ...data4ref.value]
})
watch([Xcol, Ycol, () => props.index], () => {
  data4check.value = dataCSV.value.data.map((i) => {
    return {
      name: props.allFiles[props.index]?.name,
      nm: i[XcolName.value],
      intensity: i[YcolName.value],
    }
  })
}, {
  immediate: true,
})

const updateKey = ref(0)
const figOptions: any = computed(() => ({
  type: 'line',
  data: data4show.value,
  encode: {
    x: 'nm',
    y: 'intensity',
    color: 'name',
  },
  scale: {
    x: { tickMethod: () => tickInfo(250, 650) },
    y: { nice: false, tickMethod: () => [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] },
  },
}))
</script>

<template>
  <div class="h-full w-full box-border flex justify-center items-center">
    <div class="box-border w-70% h-full flex-auto flex b-2 rd-1 b-indigo-100 b-solid">
      <g2View v-if="dataCSV.data.length > 0" :key="updateKey" :option="figOptions" />
    </div>
    <div class="flex-auto flex flex-col bg-gray-50 h-full w-30%  box-border">
      当前文件：{{ props.allFiles[props.index]?.name }}
      <NSelect v-model:value="Xcol" :options="colsOptions" size="small" :render-tag="renderTag" />
      <NSelect v-model:value="Ycol" :options="colsOptions" size="small" :render-tag="renderTag" />
      <NButton type="info" secondary size="tiny" class="flex-none text-1em m-1 ml-a " @click="updateKey++">
        重新渲染
      </NButton>
      <NButton type="info" secondary size="tiny" class="flex-none text-1em m-1 ml-a " @click="updateKey++">
        检查数据
      </NButton>
      <!-- <figCheckOption :data4check="objArr.data" /> -->
    </div>
  </div>
</template>

<style>
</style>
