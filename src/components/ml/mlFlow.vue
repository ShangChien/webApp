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
import type { dataUnimol } from '@/components/types'

const props = defineProps<{
  allFiles: { name: string; contents: string }[]
  index: number | null
}>()
const { fileType, task, models } = useMlSetting()

const currentFile = computed<{ name: string; contents: string }>(() => props.allFiles[props.index])
const allFiles = computed<{ name: string; contents: string }[]>(() => props.allFiles)

const result = ref<{ name: string; homo: number; lumo: number; eg: number; key: number }[]>([
  { key: 0, name: 'B222', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'B222', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'B222', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'B222', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
  { key: 0, name: 'E322', homo: 98, lumo: 60, eg: 70 },
])
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
    models: models.value,
    names: nameList,
    smiles: fileType.value === '*.smi' ? contentsList : null,
    molBlocks: fileType.value !== '*.smi' ? contentsList : null,
  }
})

function Inference(taskInfo: dataUnimol | Ref<dataUnimol>) {
  console.log(toValue(taskInfo))
  axios.post(
    'api/unimol',
    toValue(taskInfo),
  ).then(async (res: any) => {
    console.log(res)
  }).catch((error) => {
    console.log(error)
  })
}
</script>

<template>
  <Splitpanes class="h-full w-full text-lg default-theme " horizontal>
    <Pane min-size="20" max-size="70" size="40">
      <div class="h-full w-full flex flex-nowrap justify-start items-center box-border">
        <div class="h-full flex-auto aspect-ratio-square box-border pr-1 pb-1 bg-slate-1">
          <div v-if="currentFile?.contents" class="h-full w-full">
            <nglViewer :data="currentFile?.contents" class="w-full h-full " />
          </div>
          <div v-else class="flex flex-col justify-center items-center bg-slate-2 h-full w-full rd-1 box-border">
            <div class="i-noto-hamburger animate-bounce-alt animate-count-infinite animate-duration-1.5s text-4xl" />
            <span>no structure</span>
          </div>
        </div>
        <div class="h-full w-80% flex-auto box-border rd-1 bg-slate-1 pb-1">
          <div class="rd-1 bg-slate-2 h-full box-border flex flex-col ">
            <div class="flex flex-none justify-end items-center gap-2 p-1 pb-0 ">
              <NButton type="info" secondary class=" text-lg" @click="Inference(currentFileTaskInfo)">
                predict current
              </NButton>
              <NButton type="info" class="text-lg" @click="Inference(multiFileTaskInfo)">
                predict all
              </NButton>
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
