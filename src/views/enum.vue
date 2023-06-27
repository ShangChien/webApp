<script setup lang="ts">
import { computed, inject, onMounted, provide, ref, watch } from 'vue'
import type { Ref } from 'vue'
import {
  NAvatar, NButton, NDropdown, NGrid, NGridItem, NPageHeader,
  NSpace, NStatistic, NStep, NSteps, NSwitch, NThing,
} from 'naive-ui'
import siteSelcet from '@/components/rdkitComponent/siteSelect.vue'
import classTree from '@/components/enumMole/classTree.vue'
import combo from '@/components/enumMole/combo.vue'
import result from '@/components/enumMole/result.vue'
import type { molData } from '@/components/types'

// import localForage from "localforage";

const options = [{ label: 'Ketcher: 2.4.0' }, { label: 'RDKit: 2022.3.2' }]
const hasMounted = ref(false)
const currentStep = ref(2)// 当前步骤
const visiualBox = ref(false)// 是否显示可视化框
// 控制浮动框添加按钮的显示类型
provide('visiualBox', visiualBox)
// 位点选取组件所显示分子的id
// currentEdit.value.state:
// 0:没有操作的状态; 1:处于预览状态; 2:处于编辑状态;3:处于编辑结束
const currentEdit: Ref<{ id: number;state: number }> = inject('currentEdit')
// 注入结果
const enumResult = ref<molData[]>([])
provide('enumResult', enumResult)

const ligandsDom = ref()
const coresDom = ref()
const ligands = computed(() => hasMounted.value ? ligandsDom.value.gridData : null)
const cores = computed(() => hasMounted.value ? coresDom.value.gridData : null)
onMounted(() => {
  hasMounted.value = true
})
watch(
  [ligands, cores],
  (n, o) => {
    // console.log('old: ',o)
    // console.log('new: ',n)
  },
)
</script>

<template>
  <div>
    <NPageHeader subtitle="Molecule Edit, Substructure Search( beta 1.0 )">
      <template #title>
        <a style="text-decoration: none; color: inherit">Molecule Generate</a>
      </template>
      <template #avatar>
        <NAvatar color=" white" src="/favicon.svg" />
      </template>
      <template #extra>
        <NSpace class="mr-2">
          <NSwitch v-model:value="visiualBox">
            <template #checked>
              已显示画板
            </template>
            <template #unchecked>
              已隐藏画板
            </template>
          </NSwitch>
          <NDropdown :options="options" placement="bottom-start">
            <div class="i-noto-v1-information text-xl" />
          </NDropdown>
        </NSpace>
      </template>
      <template #default>
        <site-selcet :id="currentEdit.id" />
        <div class="b-2 rd-2 b-indigo-100 mt--3 step ">
          <div class="bg-gray-50 m-2 rd-1.5">
            <NSteps v-model:current="currentStep" class="ma-2 p-2  ">
              <NStep title="配体" description="标记配体位点和分类" />
              <NStep title="主体" description="标记主核位点和分类" />
              <NStep title="位点组合" description="位点组合" />
              <NStep title="结果展示" description="结果查看和筛选" />
            </NSteps>
          </div>
          <div class="grid grid-rows-2 items-center mr-2">
            <NButton :disabled="currentStep in [0, 1]" type="info" @click="currentStep--">
              back
            </NButton>
            <NButton :disabled="currentStep === 4" type="info" @click="currentStep++">
              {{ currentStep === 0 ? 'initialze' : 'next' }}
            </NButton>
          </div>
        </div>
      </template>
    </NPageHeader>
    <div v-show="currentStep === 1">
      <class-tree ref="ligandsDom" />
    </div>
    <div v-show="currentStep === 2">
      <class-tree ref="coresDom" />
    </div>
    <div v-show="currentStep === 3">
      <combo :cores="cores" :ligands="ligands" />
    </div>
    <div v-show="currentStep === 4">
      <result :mols="enumResult" />
    </div>
    <div v-show="![1, 2, 3, 4].includes(currentStep)">
      <NThing>
        <NGrid :cols="5">
          <NGridItem>
            <NStatistic label="主核" value="125 个" />
          </NGridItem>
          <NGridItem>
            <NStatistic label="配体" value="22 位" />
          </NGridItem>
          <NGridItem>
            <NStatistic label="结构" value="36 次" />
          </NGridItem>
          <NGridItem>
            <NStatistic label="计算次数" value="83 个" />
          </NGridItem>
        </NGrid>
      </NThing>
    </div>
  </div>
</template>

<style>
.step {
  display: grid;
  align-items: stretch;
  grid-column-gap: 0.3em;
  grid-template-columns: 1fr 80px;
}
</style>
