<script setup lang='ts'>
import { computed, ref, toValue } from 'vue'
import type { Ref } from 'vue'
import { NButton } from 'naive-ui'
import { Pane, Splitpanes } from 'splitpanes'
import axios from 'axios'
import nglViewer from '@/components/nglViewer.vue'
import mlTable from '@/components/ml/mlTable.vue'
import modelDetail from '@/components/ml/modelDetail.vue'
import { useMlSetting } from '@/stores/mlSetting'
import type { dataResults, dataUnimol } from '@/components/types'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
  index: number | null
}>()
const { fileType, task, models } = useMlSetting()

const currentFile = computed<{ name: string; contents: string }>(() => props.allFiles[props.index])
const allFiles = computed<{ name: string; contents: string }[]>(() => props.allFiles)

const result = ref<dataResults[]>([])
const currentFileTaskInfo = computed(() => ({
  atoms: null,
  coordinates: null,
  results: null,
  models: models.value,
  names: [currentFile.value.name],
  smiles: fileType.value === '*.smi' ? [currentFile.value.contents] : null,
  molBlocks: fileType.value !== '*.smi' ? [currentFile.value.contents] : null,
}))
const multiFileTaskInfo = computed(() => {
  const contentsList = allFiles.value.map(e => e.contents)
  const nameList = allFiles.value.map(e => e.name)
  return {
    atoms: null,
    coordinates: null,
    results: null,
    models: models.value,
    names: nameList,
    smiles: fileType.value === '*.smi' ? contentsList : null,
    molBlocks: fileType.value !== '*.smi' ? contentsList : null,
  }
})

const inferencing = ref(false)
function Inference(taskInfo: dataUnimol | Ref<dataUnimol>) {
  inferencing.value = true
  console.log(toValue(taskInfo))
  axios.post(
    'api/unimol',
    toValue(taskInfo),
  ).then(async (res: any) => {
    result.value = res2Obj(res.data)
    inferencing.value = false
  }).catch((error) => {
    console.log(error)
    inferencing.value = false
  })
}

function res2Obj(res: dataUnimol) {
  const nameList = res.names
  const keys = Object.keys(res.results)
  const outObj = nameList.map((el, index) => {
    const newObj = keys.reduce((obj, key) => {
      const val = res.results[key][index]
      obj[key] = Array.isArray(val) ? val : 27.21138 * (val as number)
      return obj
    }, {})
    return {
      name: el.endsWith('.mol') ? el.slice(0, -4) : el,
      key: index,
      ...newObj,
    }
  })
  return outObj
}
</script>

<template>
  <Splitpanes class="h-full w-full text-lg default-theme " horizontal>
    <Pane min-size="20" max-size="70" size="40">
      <div class="h-full w-full flex flex-nowrap justify-start items-center box-border">
        <div class="h-full flex-auto aspect-ratio-square box-border pr-1 pb-1 bg-slate-1 flex">
          <div v-if="currentFile?.contents" class="h-full w-full">
            <nglViewer :data="currentFile?.contents" class="w-full h-full " />
          </div>
          <div v-else class="ma text-slate text-center">
            <span>3D Viewer</span>
            <div class="i-carbon-face-dizzy text-4xl ma" />
            <span>no data</span>
          </div>
        </div>
        <div class="h-full w-80% flex-auto box-border rd-1 bg-slate-1 pb-1">
          <div class="rd-1 bg-slate-2 h-full box-border flex flex-col ">
            <div class="flex-none h-38px">
              <div v-if="!inferencing" class="flex justify-end items-center gap-2 p-1 pb-0 ">
                <NButton type="info" secondary class=" text-lg" @click="Inference(currentFileTaskInfo)">
                  predict current
                </NButton>
                <NButton type="info" class="text-lg" @click="Inference(multiFileTaskInfo)">
                  predict all
                </NButton>
              </div>
              <div v-else class=" h-full flex justify-end items-center gap-2 p-1 pb-0 pr-2 text-violet-5">
                <div class="i-svg-spinners-bars-scale" />
                inferencing...
              </div>
            </div>
            <div class="p-1 flex-auto h-80% box-border">
              <div class="rd-1 bg-slate-1 p-2 flex-auto h-full box-border">
                <modelDetail v-model:task="task" v-model:fileType="fileType" v-model:models="models" class="box-border" />
              </div>
            </div>
          </div>
        </div>
      </div>
    </Pane>
    <Pane max-size="70">
      <mlTable class="flex-auto bg-slate-1 py-1 h-full" :data="result" />
    </Pane>
  </Splitpanes>
</template>

<style>
@import 'splitpanes/dist/splitpanes.css';

.splitpanes__pane {
  display: flex;
  justify-content: center;
  align-items: center;
}
</style>
