<script setup lang="ts">
import { ref, onMounted,provide,computed } from "vue";
import { NStep,NSteps,NCard,NButton,NButtonGroup,NIcon,NSpace,NDivider,NSwitch,NPageHeader,
NGrid,NGridItem,NStatistic,NAvatar,NDropdown,NThing } from "naive-ui"
import siteSelcet from "@/components/rdkitComponent/siteSelect.vue"
import classTree from "@/components/enumMole/classTree.vue"
import type { molData } from "@/components/types"
import type { ComputedRef } from "vue";
import { useEnumStore } from '@/stores/enumStore'
//import localForage from "localforage";
//relavent lib
const options = [{ label: "Ketcher: 2.4.0" }, { label: "RDKit: 2022.3.2" }];
const currentStep=ref(0)//当前步骤
const visiualBox=ref(true)//是否显示可视化框
provide('visiualBox',visiualBox)
//控制浮动框添加按钮的显示类型
const molAdd = computed(()=>{
  if (currentStep.value==1) {
    return {molType:'ligand',info:'添加配体'}
  } else if (currentStep.value==2) {
    return {molType:'core',info:'添加主核'}
  } else {
    return {molType:'mole',info:'添加分子'}
  }
})
provide('molAdd',molAdd)

const cores=ref<molData[]|any>()
const ligands=ref<molData[]|any>()

const enumStore = useEnumStore()
//监听pinia状态action的调用，数据持久化保存到本地

onMounted(()=>{

})

</script>

<template>
<div >
<n-page-header subtitle="Molecule Edit, Substructure Search( beta 1.0 )">
  <template #title>
    <a style="text-decoration: none; color: inherit">Molecule Generate</a>
  </template>
  <template #avatar>
    <n-avatar
      color=" white"
      src="/favicon.svg"
    />
  </template>
  <template #extra>
    <n-space class="mr-2">
      <n-switch v-model:value="visiualBox" >
        <template #checked>
          已显示画板
        </template>
        <template #unchecked>
          已隐藏画板
        </template>
      </n-switch>
      <n-dropdown :options="options" placement="bottom-start">
        <div class="i-noto-v1-information text-xl"></div>
      </n-dropdown>
    </n-space>
  </template>
  <template #default>
    <site-selcet />
    <n-card embedded class="mt--3">
      <n-steps v-model:current="currentStep">
        <n-step title="配体" description="标记配体位点和分类" />
        <n-step title="主体" description="标记主核位点和分类" />
        <n-step title="位点组合" description="位点组合" />
    		<n-step title="结果展示" description="结果查看和筛选" />
      </n-steps>
    	<n-space style="padding-top:10px" >
        <n-button @click="currentStep--" :disabled="currentStep in [0,1]" type="info">
          back
        </n-button>
    		<n-button @click="currentStep++" :disabled="currentStep===4" type="info">
          {{currentStep === 0 ? 'initialze' : 'next'}}
        </n-button>
    	</n-space>
    </n-card>
  </template>  
</n-page-header>

<div v-if="currentStep === 1" style="padding-top:10px;">
  <class-tree style="width:20%"/>
  <!-- <grid-page :mollist="ligands" /> -->
</div>
<div v-else-if="currentStep === 2">
  主核
  <!-- <grid-page :mollist="cores" /> -->
</div>
<div v-else-if="currentStep === 3">
  组合设置
</div>
<div v-else-if="currentStep === 4">
  结果页
</div>
<div v-else>
<n-thing>
  <n-grid :cols="5" >
    <n-grid-item>
      <n-statistic label="主核" value="125 个" />
    </n-grid-item>
    <n-grid-item>
      <n-statistic label="配体" value="22 位" />
    </n-grid-item>
    <n-grid-item>
      <n-statistic label="结构" value="36 次" />
    </n-grid-item>
    <n-grid-item>
      <n-statistic label="计算次数" value="83 个" />
    </n-grid-item>
  </n-grid>
</n-thing>

</div>
</div>
</template>