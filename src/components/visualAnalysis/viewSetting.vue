<script setup lang='ts'>
import { computed, h, reactive, ref, watch } from 'vue'
import { NButton, NCollapse, NCollapseItem, NScrollbar, NSelect, NSlider } from 'naive-ui'
import type { ScaleItem, ViewType, itemOption } from './types'
import { coordOption, fig, scaleAxisType, scaleType, view } from './types'
import { toOptionsArray } from './utils'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
}>()
const emit = defineEmits<{ updateChart: [data: ViewType] }>()
const fileIndex = defineModel<number | null>('index', { default: 0 })
const config = defineModel<ViewType>('config')

const fileOptions = computed(() => props.allFiles.map((el, index) => ({ label: el.name, value: index })))
const columns = computed(() => {
  const file = props.allFiles[fileIndex.value]
  let cols
  if (file === null || file === undefined || Object.keys(file).length === 0) {
    return []
  } else {
    if (file.name.endsWith('csv') || file.name.endsWith('CSV')) {
      cols = file.contents.split('\n')[0].split(',')
    } else if (file.name.endsWith('json') || file.name.endsWith('JSON')) {
      cols = Object.keys(JSON.parse(file.contents)[0])
    }
    return cols
  }
  // return ['height', 'width', 'length']
})

const viewView = reactive({
  value: view.view,
  options: toOptionsArray(view),
})
const viewCoord = reactive({
  value: coordOption.cartesian,
  options: toOptionsArray(coordOption),
})
const scaleConfig: any = Object.entries(scaleAxisType).reduce((obj, [key, _value]) => {
  obj[key] = {
    type: {
      value: scaleType.linear,
      options: toOptionsArray(scaleType),
    },
    range: {
      value: [0, 100],
      options: null,
    },
    domain: {
      value: [0, 100],
      options: null,
    },
    nice: {
      value: 'True',
      options: ['True', 'False'].map(x => ({ value: x, label: x })),
    },
    padding: {
      value: 0,
      options: [5, 1, 0, -1, -5].map(x => ({ value: x, label: x })),
    },
  }
  return obj
}, reactive({}))

const XYs = ref<any[]>([{
  type: {
    value: fig.point,
    options: toOptionsArray(fig),
  },
  encode: {
    x: {
      value: columns.value[0],
      options: columns.value.map(x => ({ value: x, label: x })),
    },
    y: {
      value: columns.value[1],
      options: columns.value.map(x => ({ value: x, label: x })),
    },
    color: {
      value: columns.value[1],
      options: columns.value.map(x => ({ value: x, label: x })),
    },
  },
  scale: scaleConfig,
}])

watch(XYs, () => {
  console.log(scaleConfig)
}, { deep: true })

function logic2component(XY: any) {
  return h('div', { class: 'grid gap-1 box-border items-center', style: { 'grid-template-columns': '60px auto' } },
    Object.entries(XY).map((item: any) => {
      if ('options' in item[1]) {
        return [
          h('div', { class: 'bg-slate-1 rd-1 h-full w-full items-center justify-center flex' }, item[0]),
          Array.isArray(item[1].options)
            ? h(NSelect, {
              'size': 'small',
              'value': item[1].value,
              'options': item[1].options,
              'onUpdate:value': (value: any) => { item[1].value = value },
            }, null)
            : h(NSlider, {
              'range': true,
              'step': 1,
              'value': item[1].value,
              'onUpdate:value': (value: any) => { item[1].value = value },
            }, null),
        ]
      } else {
        return [h('div', { class: 'bg-slate-1 rd-1 h-full w-full items-center justify-center flex' }, item[0]), logic2component(item[1])]
      }
    }).flat(),
  )
}
</script>

<template>
  <div class="flex-auto box-border h-full flex flex-col">
    <NButton type="info" secondary size="tiny" class="flex-none text-1em m-1 ml-a " @click="emit('updateChart', config)">
      重新渲染
    </NButton>
    <NScrollbar class="max-h-full pr-1 box-border flex flex-auto">
      <NCollapse>
        <NCollapseItem title="视图布局" name="1">
          <div class="grid grid-cols-2 gap-1 box-border m-2" style="grid-template-columns:100px auto">
            <div class="bg-slate-1  rd-1 h-full w-full items-center justify-center flex">
              数据源：
            </div>
            <NSelect
              v-model:value="fileIndex"
              size="small"
              :options="fileOptions"
              :disabled="fileOptions.length < 1"
              :placeholder="(fileOptions.length < 1) ? 'no data' : props.allFiles[fileIndex].name"
            />
            <div class="bg-slate-1  rd-1 h-full w-full items-center justify-center flex">
              布局类型：
            </div>
            <NSelect v-model:value="viewView.value" size="small" :options="viewView.options" />
            <div class="bg-slate-1  rd-1 h-full w-full items-center justify-center flex">
              坐标系类型：
            </div>
            <NSelect v-model:value="viewCoord.value" size="small" :options="viewCoord.options" />
          </div>
        </NCollapseItem>
        <NCollapseItem title="子图设置" name="2">
          <div v-for="(xy, key) in XYs" :key="key" class="m-2">
            <component :is="logic2component(xy)" />
          </div>
        </NCollapseItem>
      </NCollapse>
    </NScrollbar>
  </div>
</template>

<style>
</style>
