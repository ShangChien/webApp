<script setup lang='ts'>
import { computed, reactive, ref } from 'vue'
import { NButton, NCollapse, NCollapseItem, NScrollbar, NSelect } from 'naive-ui'
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
      value: null,
      option: {
        min: 0,
        max: 100,
      },
    },
    domain: {
      value: null,
      option: {
        min: 0,
        max: 100,
      },
    },
    nice: {
      value: true,
      option: [true, false],
    },
    padding: {
      value: 0.01,
      option: [true, false],
    },
  }
  return obj
}, reactive({}))

const XYs = ref<itemOption[]>([{
  type: fig.point,
  encode: {
    x: 'weight',
    y: 'height',
    color: 'gender',
  },
  scale: scaleConfig,
}])

const localConfig = computed<ViewType>(() => ({
  type: viewView.value,
  data: JSON.parse(props.allFiles[fileIndex.value].contents),
  coordinate: { type: viewCoord.value },
  scales: {
    x: { nice: true },
    y: { nice: true },
  },
  children: XYs.value,
}))
</script>

<template>
  <div>
    <NButton type="info" secondary size="tiny" class=" text-1em float-right m-1" @click="emit('updateChart', config)">
      重新渲染
    </NButton>
    <NScrollbar class="max-h-full pr-1 box-border">
      <NCollapse class="m-1">
        <NCollapseItem title="视图布局" name="1">
          <div class="grid grid-cols-2 gap-1 box-border m-1">
            <div>数据源：</div>
            <NSelect
              v-model:value="fileIndex"
              size="small"
              :options="fileOptions"
              :disabled="fileOptions.length < 1"
              :placeholder="(fileOptions.length < 1) ? 'no data' : props.allFiles[fileIndex].name"
            />
            <div>布局类型：</div>
            <NSelect v-model:value="viewView.value" size="small" :options="viewView.options" />
            <div>坐标系类型：</div>
            <NSelect v-model:value="viewCoord.value" size="small" :options="viewCoord.options" />
          </div>
          <NCollapse>
            <NCollapseItem title="坐标轴" name="1-1">
              <NCollapse>
                <NCollapseItem title="x轴" name="1-1">
                  <div class="grid grid-cols-2 gap-1 box-border m-1">
                    <div>类型：</div>
                    <NSelect
                      v-model:value="fileIndex"
                      size="small"
                      :options="fileOptions"
                    />
                    <div>布局类型：</div>
                    <NSelect v-model:value="viewLayout" size="small" :options="viewLayoutOptions" />
                    <div>坐标系类型：</div>
                    <NSelect v-model:value="coordLayout" size="small" :options="coordLayoutOptions" />
                  </div>
                </NCollapseItem>
                <NCollapseItem title="y轴" name="1-2">
                  <div>快速通过</div>
                </NCollapseItem>
              </NCollapse>
            </NCollapseItem>
          </NCollapse>
        </NCollapseItem>
        <NCollapseItem title="子图设置" name="2">
          <div>布局类型：</div>
          <NSelect v-model:value="viewLayout" size="small" :options="viewLayoutOptions" />
          <div>布局类型：</div>
          <NSelect v-model:value="viewLayout" size="small" :options="viewLayoutOptions" />
        </NCollapseItem>
      </NCollapse>
    </NScrollbar>
  </div>
</template>

<style>
</style>
