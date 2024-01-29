<script setup lang='ts'>
import { computed, defineAsyncComponent, ref, watch } from 'vue'
import { NButton } from 'naive-ui'
import type { spectrumData } from './dataCheckStore'
import { file2Data, normalize, tickInfo } from '@/components/visualAnalysis/utils'
import figCheckOption from '@/components/dataCheck/figCheckOption.vue'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
  index: number | null
}>()

const g2View = defineAsyncComponent(() => import('@/components/visualAnalysis/g2View.vue'))

const objArr = computed(() => {
  const _out = file2Data(props.allFiles[props.index])
  _out.data = normalize(_out.data)
  _out.data = _out.data.map(item => ({
    nm: item[_out.cols[0]],
    intensity: item[_out.cols[1]],
    name: props.allFiles[props.index].name.split('.')[0],
  }))
  _out.cols = Object.keys(_out.data[0] ?? {})
  return _out
})

const figOptions: any = computed(() => ({
  type: 'line',
  data: objArr.value.data,
  encode: {
    x: objArr.value.cols[0],
    y: objArr.value.cols[1],
    color: 'name',
  },
  scale: {
    x: { tickMethod: () => tickInfo(250, 650) },
    y: { nice: false, tickMethod: () => [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] },
  },
}))
const updateKey = ref(0)

watch(figOptions, () => {
  console.log(props)
}, {
  deep: true,
  immediate: true,
})
</script>

<template>
  <div class="h-full w-full box-border flex justify-center items-center">
    <div class="box-border w-70% h-full flex-auto flex b-2 rd-1 b-indigo-100 b-solid">
      <g2View v-if="objArr.data.length > 0" :key="updateKey" :option="figOptions" />
    </div>
    <div class="flex-auto flex flex-col bg-gray-50 h-full w-30%  box-border">
      当前文件：{{ props.allFiles[props.index]?.name }}
      <NButton type="info" secondary size="tiny" class="flex-none text-1em m-1 ml-a " @click="updateKey++">
        重新渲染
      </NButton>
      <NButton type="info" secondary size="tiny" class="flex-none text-1em m-1 ml-a " @click="updateKey++">
        检查数据
      </NButton>
      <figCheckOption :data4check="objArr.data" />
    </div>
  </div>
</template>

<style>
</style>
