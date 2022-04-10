<script setup lang="ts">
import { NSpace,NPageHeader,NGrid,NGridItem,NStatistic,NBreadcrumb,NMessageProvider,
NBreadcrumbItem,NButton,NDropdown,NAvatar,useMessage,NIcon,NSwitch,NCollapse,NCard,
NCollapseItem,NInputGroup,NInput,NDivider} from 'naive-ui'
import { ColorPaletteOutline } from '@vicons/ionicons5'
import initKetcher from '../components/initKetcher.vue';
import editableRdkit from '@/components/editableRdkit.vue'
import svgRdkit from '@/components/svgRdkit.vue'
import cardRdkit from '@/components/cardRdkit.vue'
import { reactive,onMounted,ref, nextTick, watch } from 'vue';
import type { molData } from '@/components/types'

const inputText=ref('')
const showAbout=ref(false)
const editRdkitKey=ref(1)
const active = ref(true)
const options=[{label: 'Ketcher: 2.4.0'},{label: 'RDKit: 2021.9.5'}]

function drawMol(){
  initMol.smiles=inputText.value
  //console.log(initMol)
}

const initSmile=ref<any>('CC(=O)Oc1ccccc1C(=O)O')
const initMol:molData=reactive({
    smiles:initSmile.value,
    qsmiles:'*~*',
    width:500,
    height:500,
    addAtomIndices:true,
})

const tmpSmile=reactive({
  smiles:initMol.smiles as any,
  atoms:[] as any,
  bonds:[] as any,
})
function acceptMol(mol:molData){
  tmpSmile.smiles=mol.smiles
  tmpSmile.atoms=mol.atoms
  tmpSmile.bonds=mol.bonds
  console.log(tmpSmile)
}

const mol4enum=reactive({
  core:[] as Array<molData>,
  ligand:[] as Array<molData>
})
//addcore
function addCore(){
  editRdkitKey.value++
  const mol={
    smiles:'',
    atoms:[],
    bonds:[],
  }
  //判断存在形式
  if((tmpSmile.smiles!='')&&((JSON.stringify(mol4enum.core).indexOf(JSON.stringify(tmpSmile))==-1))){
    mol.smiles=tmpSmile.smiles
    mol.atoms=tmpSmile.atoms
    mol.bonds=tmpSmile.bonds
    mol4enum.core.push(JSON.parse(JSON.stringify(mol)))
    //tmpSmile.smiles=''
    //tmpSmile.atoms=[]
    //tmpSmile.bonds=[]
  }
  console.log(mol4enum.core)
}
//addLigand
function addLigand(){
  editRdkitKey.value++
  const mol={
    smiles:'',
    atoms:[],
    bonds:[],
  }
  //判断存在形式
  if((tmpSmile.smiles!='')&&((JSON.stringify(mol4enum.ligand).indexOf(JSON.stringify(tmpSmile.smiles))==-1)
  ||(JSON.stringify(mol4enum.ligand).indexOf(JSON.stringify(tmpSmile.bonds))==-1)
  ||(JSON.stringify(mol4enum.ligand).indexOf(JSON.stringify(tmpSmile.atoms))==-1))){
    mol.smiles=tmpSmile.smiles
    mol.atoms=tmpSmile.atoms
    mol.bonds=tmpSmile.bonds
    mol4enum.ligand.push(JSON.parse(JSON.stringify(mol)))
    //tmpSmile.smiles=''
    //tmpSmile.atoms=[]
    //tmpSmile.bonds=[]
  }
  console.log(mol4enum.ligand)
}


onMounted(()=>{
})

watch(
  ()=>initMol.smiles,
  (val)=>{
  tmpSmile.smiles=val
})
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
      <n-avatar @click='showAbout=!showAbout'
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
                <n-input v-model:value="inputText"
                         style="font-size:20px"
                         type="text"
                         size="large"
                         clearable
                         placeholder="type smiles here"/>
                <n-button size="large" style="font-size:20px" color="#c471ed" @click="drawMol">
                          <template #icon>
                            <n-icon>
                              <ColorPaletteOutline />
                            </n-icon>
                          </template>
                   绘制
                </n-button>
            </n-input-group>
            <n-space justify="space-around">
              <!--rdkit-sub v-bind="initMol"/-->
              <div style="width:500px ;height:500px">
              <editable-rdkit :key="editRdkitKey" v-bind="initMol" @update-mol="acceptMol"/>
              </div>
              <n-space align="stretch" :size="60">
                <n-button size="medium" strong secondary round type="info" @click="addCore" >添加主核</n-button>
                <n-button size="medium" strong secondary round type="info" @click="addLigand" >添加配体</n-button>
                <n-button szie="large" round color="#ff69b4" >开始枚举</n-button>
              </n-space>  
            </n-space>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="Ketcher画板" display-directive="show" name="1">
          <init-ketcher class="collapse" @update-smiles="(smiles)=> inputText=smiles"></init-ketcher>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="主核结构" display-directive="show" name="1">
          <n-grid :cols="4" x-gap="8" y-gap="8">
            <n-grid-item
              v-for="item in mol4enum.core"
              :key="mol4enum.core.indexOf(item)">
              <card-rdkit v-bind="item" style="width:100%;"/>
            </n-grid-item>
          </n-grid>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="配体结构" display-directive="show" name="1">
          <n-grid :cols="6" x-gap="8" y-gap="8">
            <n-grid-item
              v-for="item in mol4enum.ligand"
              :key="mol4enum.ligand.indexOf(item)">
              <card-rdkit v-bind="item"/>
            </n-grid-item>
          </n-grid>
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
