<script setup lang="ts">
import { ref, onMounted,provide,computed,watch } from "vue";
import { NStep,NSteps,NButton,NSpace,NSwitch,NPageHeader,
NGrid,NGridItem,NStatistic,NAvatar,NDropdown,NThing } from "naive-ui"
import siteSelcet from "@/components/rdkitComponent/siteSelect.vue"
import classTree from "@/components/enumMole/classTree.vue"
import combo from "@/components/enumMole/combo.vue"
import result from "@/components/enumMole/result.vue"
import gridPage from "@/components/rdkitComponent/gridPage.vue";
import type { molData } from "@/components/types"
import { useEnumStore } from '@/stores/enumStore'

//import localForage from "localforage";

const options = [{ label: "Ketcher: 2.4.0" }, { label: "RDKit: 2022.3.2" }];
const hasMounted=ref(false)
const currentStep=ref(2)//当前步骤
const visiualBox=ref(false)//是否显示可视化框
//控制浮动框添加按钮的显示类型
provide('visiualBox',visiualBox)
//位点选取组件所显示分子的id
//currentEdit.value.state:
//0:没有操作的状态; 1:处于预览状态; 2:处于编辑状态;3:处于编辑结束
const currentEdit = ref<{id:number;state:number}>({id:0,state:0})
provide('currentEdit',currentEdit)
//注入结果
const enumResult=ref<molData[]>([])
provide('enumResult',enumResult)

const ligandsDom=ref()
const coresDom=ref()
const ligands=computed(()=>hasMounted.value?ligandsDom.value.gridData:null)
const cores=computed(()=>hasMounted.value?coresDom.value.gridData:null)
onMounted(()=>{
  hasMounted.value=true
})
watch(
  [ligands,cores],
  (n,o)=>{
    //console.log('old: ',o)
    //console.log('new: ',n)
  }
)

</script>

<template>
<div >
<n-page-header subtitle="Molecule Edit, Substructure Search( beta 1.0 )">
  <template #title>
    <a style="text-decoration: none; color: inherit">Molecule Generate</a>
  </template>
  <template #avatar>
    <n-avatar color=" white" src="/favicon.svg" />
  </template>
  <template #extra>
    <n-space class="mr-2">
      <n-switch v-model:value="visiualBox" >
        <template #checked>已显示画板</template>
        <template #unchecked> 已隐藏画板</template>
      </n-switch>
      <n-dropdown :options="options" placement="bottom-start">
        <div class="i-noto-v1-information text-xl"></div>
      </n-dropdown>
    </n-space>
  </template>
  <template #default>
    <site-selcet :id='currentEdit.id' />
    <div class="b-2 rd-2 b-indigo-100 mt--3 step ">
      <div class="bg-gray-50 m-2 rd-1.5" ><n-steps v-model:current="currentStep" class="ma-2 p-2  ">
        <n-step title="配体" description="标记配体位点和分类" />
        <n-step title="主体" description="标记主核位点和分类" />
        <n-step title="位点组合" description="位点组合" />
    		<n-step title="结果展示" description="结果查看和筛选" />
      </n-steps></div>
    	<div class="grid grid-rows-2 items-center mr-2" >
        <n-button  @click="currentStep--" :disabled="currentStep in [0,1]" type="info">
          back
        </n-button>
    		<n-button @click="currentStep++" :disabled="currentStep===4" type="info">
          {{currentStep === 0 ? 'initialze' : 'next'}}
        </n-button>
    	</div>
    </div>
  </template>  
</n-page-header>
<div v-show="currentStep === 1" >
  <class-tree  ref='ligandsDom'/>
</div>
<div v-show="currentStep === 2">
  <class-tree  ref='coresDom'/>
</div>
<div v-show="currentStep === 3">
  <combo :cores="cores" :ligands="ligands"></combo> 
</div>
<div v-show="currentStep === 4">
  <result :mols="enumResult" />
</div>
<div v-show="![1,2,3,4].includes(currentStep)">
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
<style>
.step {
  display: grid;
  align-items: stretch;
  grid-column-gap: 0.3em;
  grid-template-columns: 1fr 80px;
}
</style>