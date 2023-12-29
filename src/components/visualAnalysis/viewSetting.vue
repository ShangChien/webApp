<script setup lang='ts'>
import { computed, defineAsyncComponent, h, reactive, ref, watch } from 'vue'
import { NButton, NCollapse, NCollapseItem, NScrollbar, NSelect, NSlider } from 'naive-ui'
import type { ViewType, itemOption } from './types'
import { coordOption, fig, scaleAxisType, scaleType, view } from './types'
import { readCsv, toOptionsArray } from './utils'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
}>()
const fileIndex = defineModel<number>('index')

const g2View = defineAsyncComponent(() => import('./g2View.vue'))
const fileOptions = computed(() => props.allFiles.map((el, index) => ({ label: el.name, value: index })))
const dataObj = computed(() => {
  const file = props.allFiles[fileIndex.value]
  if (file === null || file === undefined || Object.keys(file).length === 0) {
    return { cols: [], data: [] }
  } else {
    let cols: string[]
    let data: object[]
    if (file.name.endsWith('csv') || file.name.endsWith('CSV')) {
      data = readCsv(file.contents)
      cols = file.contents.split('\n')[0].split(',')
    } else if (file.name.endsWith('json') || file.name.endsWith('JSON')) {
      data = JSON.parse(file.contents)
      cols = Object.keys(data[0])
    }
    return { cols, data }
  }
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
}, {})

function createEncode(cols: any[]) {
  const options = cols.map((x: any) => ({ value: x, label: x }))
  return {
    x: { value: cols[0], options },
    y: { value: cols[1], options },
    color: { value: cols[1], options },
  }
}
const XYs = ref<any[]>([{
  type: {
    value: fig.point,
    options: toOptionsArray(fig),
  },
  encode: createEncode(dataObj.value.cols),
  scale: scaleConfig,
}])
watch(dataObj, (val) => {
  XYs.value.forEach((el) => {
    el.encode = createEncode(val.cols)
  })
})

const options = computed<ViewType>(() => {
  const children = XYs.value.map((el: any) => {
    const obj: itemOption = {
      type: el.type.value,
      encode: {
        x: el.encode.x.value,
        y: el.encode.y.value,
        color: el.encode.color.value,
      },
      // scale: Object.entries(el.scale).reduce((obj1, [key, value]: any) => {
      //   obj1[key] = Object.keys(value).reduce((obj2, k) => {
      //     obj2[k] = value[k].value
      //     return obj2
      //   }, {})
      //   return obj1
      // }, {}),
    }
    return obj
  })
  return {
    type: viewView.value,
    data: dataObj.value.data,
    coordinate: { type: viewCoord.value },
    children,
  }
})
const updateKey = ref(0)

watch(options, () => {
  console.log(options.value)
})

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
  <div class="h-full w-full box-border flex justify-center items-center">
    <div class="box-border w-70% h-full flex-auto flex b-2 rd-1 b-indigo-100 b-solid">
      <g2View :key="updateKey" :option="options" />
    </div>
    <div class="flex-auto flex flex-col bg-gray-50 h-full w-30%  box-border">
      <NButton type="info" secondary size="tiny" class="flex-none text-1em m-1 ml-a " @click="updateKey++">
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
              <component :is="logic2component(xy)" :key="XYs" />
            </div>
          </NCollapseItem>
        </NCollapse>
      </NScrollbar>
    </div>
  </div>
</template>

<style>
</style>
