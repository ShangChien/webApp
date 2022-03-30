<script setup lang="ts">
import { NSpace,NPageHeader,NGrid,NGridItem,NStatistic,NBreadcrumb,NMessageProvider,
NBreadcrumbItem,NButton,NDropdown,NAvatar,useMessage,NIcon,NSwitch,NCollapse,NCard,
NCollapseItem,NInputGroup,NInput,NDivider} from 'naive-ui'
import { ColorPaletteOutline } from '@vicons/ionicons5'
import initKetcher from '../components/initKetcher.vue';
import rdkitSub from '@/components/rdkitSub.vue'
import { reactive,onMounted,ref } from 'vue';
import type { molData } from '@/components/types'

const initMol:molData=reactive({
    smiles:'CC(=O)Oc1ccccc1C(=O)O',
    qsmiles:'c1ccccc1',
    width:500,
    height:500,
    addAtomIndices:true,
    addBondIndices:true,
})
const initSmile=ref<any>(null)
const showAbout=ref(false)
const active = ref(true)
const options=[{label: 'Ketcher: 2.4.0'},{label: 'RDKit: 2021.9.5'}]
const message = useMessage()
function createMessage(){
  message.warning("I never needed anybody's help in any way")
}
function toggleShowAbout(){
  showAbout.value=!showAbout.value
}
function updateSmiles(smiles:string){
  initSmile.value=smiles
  console.log(initSmile.value)
}
function testsvg(){
  console.log('test svg')
}

</script>
<template>
  <n-page-header subtitle="分子结构枚举，子结构筛选 ( beta 1.0 )">
    <n-grid :cols="5" v-if="showAbout">
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
    <n-divider />
    <template #title>
      <n-space>
        <a style="text-decoration: none; color: inherit"
          >分子生成器</a>
        <n-switch v-model:value="active" size="medium"  />
      </n-space>
    </template>
    <template #avatar>
      <n-avatar @click='toggleShowAbout'
       color=" white" src="/picture/enum.svg"/>
    </template>
    <template #extra>
      <n-space>
        <n-dropdown :options="options" placement="bottom-start">
          <n-button text>· · ·</n-button>
        </n-dropdown>
      </n-space>
    </template>
  </n-page-header>
  <n-grid x-gap="40" y-gap="20" :cols="2" v-show="active" >
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title=" RDKit位点标注" display-directive="show"  class="collapse" name="1">
            <n-input-group class="inputG" >
                <n-input v-model:value="initSmile"
                         style="font-size:20px"
                         type="text"
                         size="large"
                         clearable
                         placeholder="type smiles here" 
                         />
                <n-button size="large" style="font-size:20px" color="#c471ed">
                          <template #icon>
                            <n-icon>
                              <ColorPaletteOutline />
                            </n-icon>
                          </template>
                   绘制
                </n-button>
            </n-input-group>
            <n-space justify="space-around">
              <rdkit-sub v-bind="initMol"></rdkit-sub>
              <n-space align="stretch" vertical :size="60">
                <br>
                <n-button size="medium" strong secondary round type="info">添加主核</n-button>
                <n-button size="medium" strong secondary round type="info">添加配体</n-button>
                <n-button szie="large" round color="#ff69b4" >开始枚举</n-button>
              </n-space>  
            </n-space>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="Ketcher画板" display-directive="show" name="1">
          <init-ketcher class="collapse" @update-smiles="(smiles)=> initSmile=smiles"></init-ketcher>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="主核结构" display-directive="show" name="1">
          <svg viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">
            <circle cx="100" cy="100" r="100" @click="testsvg" />
          </svg>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="配体结构" display-directive="show" name="1">
          none
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
  </n-grid>
  
</template>
<style>
.collapse{
  padding: 5px;
  width:100%;
}
.inputG{
  width:100%;
  position: relative;
  border-radius: 3px 3px 3px 3px;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0px 1px 1px 0px rgba(0, 0, 0, 0.1);
}
</style>
background-image: linear-gradient(25deg, #d24d53, #b78563, #8daf75, #27d587)