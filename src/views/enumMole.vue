<script setup lang="ts">
import {
  NAvatar,
  NButton,
  NCard,
  NCollapse,
  NCollapseItem,
  NDivider,
  NDropdown,
  NEmpty,
  NGradientText,
  NGrid,
  NGridItem,
  NIcon,
  NInput,
  NInputGroup,
  NModal,
  NPageHeader,
  NScrollbar,
  NSpace,
  NSpin,
  NStatistic,
  NSwitch,
  NThing,
} from 'naive-ui'
import { useClipboard } from '@vueuse/core'

// import editableRdkit from "@/components/editableRdkit.vue";
import {
  computed, defineAsyncComponent, reactive, ref,
  watch,
} from 'vue'
import axios from 'axios'
import Grid from 'vue-virtual-scroll-grid'
import type { molData } from '@/components/types'
import initKetcher from '@/components/ketcher/initKetcher.vue'

const cardRdkit = defineAsyncComponent({
  loader: () => import('@/components/rdkitComponent/cardRdkit.vue'),
  loadingComponent: NSpin,
  delay: 2000, // +600*Math.random(),
  errorComponent: NEmpty,
  timeout: 3000,
  suspensible: false,
})
const editableRdkit = defineAsyncComponent({
  loader: () => import('@/components/rdkitComponent/editableRdkit.vue'),
  loadingComponent: NSpin,
  delay: 2000, // +600*Math.random(),
  errorComponent: NEmpty,
  timeout: 3000,
  suspensible: false,
})

const resultData: any = ref()

// virtualGrid
const VLength = computed(() => resultData.value?.data.message.length)
const VPageSize = ref<number>(20)
const pageProvider: any = (pageNumber: number, pageSize: number) => (
  new Promise((resolve) => {
    if (resultData.value?.data.message.length < 1) {
      setTimeout(
        () => resolve(
          new Array(pageSize).fill('Loaded Item'),
        ),
        100,
      )
    } else {
      setTimeout(
        () => resolve(
          resultData.value?.data.message.slice(
            pageNumber * pageSize, (pageNumber + 1) * pageSize,
          ),
        ),
        100,
      )
    }
  })
)

const inputText = ref('')
const showAbout = ref(false)
const showModal = ref(false)
const editRdkitKey = ref(1)
const showLigandCore = ref(true)
const options = [{ label: 'Ketcher: 2.4.0' }, { label: 'RDKit: 2022.3.2' }]
const initMol: molData = reactive({
  smiles: 'CC(=O)Oc1ccccc1C(=O)O',
  qsmiles: '*~*',
  atoms: {},
  bonds: {},
})

const { copy } = useClipboard()

const tmpSmile = reactive({
  smiles: initMol.smiles,
  atoms: {} as { [key: number]: number[] } | undefined,
  bonds: {} as { [key: number]: number[] } | undefined,
})

function drawMol() {
  initMol.smiles = inputText.value
  tmpSmile.atoms = {}
  tmpSmile.bonds = {}
  // console.log(initMol)
}

function acceptMol(mol: molData) {
  tmpSmile.smiles = mol.smiles
  tmpSmile.atoms = mol.atoms
  tmpSmile.bonds = mol.bonds
  console.log(tmpSmile)
}

const mol4enum = reactive({
  core: [] as Array<molData>,
  ligand: [] as Array<molData>,
})

function enumMol() {
  axios.post('/api/enum', mol4enum)
    .then((response: any) => {
      resultData.value = response
      console.log(response)
    })
    .catch((error) => {
      console.log(error)
    })
}

// addcore
let id_core = 1
function addCore() {
  // editRdkitKey.value++ 是否刷新画板
  const mol: molData = {
    id: id_core++,
    smiles: '',
    atoms: [],
    bonds: [],
  }
  // 判断存在形式
  if (
    tmpSmile.smiles !== ''
    && (!JSON.stringify(mol4enum.core).includes(JSON.stringify(tmpSmile.smiles))
      || !JSON.stringify(mol4enum.core).includes(`"bonds":${JSON.stringify(tmpSmile.bonds)}`)
      || !JSON.stringify(mol4enum.core).includes(`"atoms":${JSON.stringify(tmpSmile.atoms)}`))
  ) {
    mol.smiles = tmpSmile.smiles
    mol.atoms = tmpSmile.atoms
    mol.bonds = tmpSmile.bonds
    mol4enum.core.push(JSON.parse(JSON.stringify(mol)))
    // tmpSmile.smiles=''
    // tmpSmile.atoms=[]
    // tmpSmile.bonds=[]
  }
  console.log(mol4enum.core)
}
// addLigand
let id_ligand = 1
function addLigand() {
  // editRdkitKey.value++
  const mol: molData = {
    id: id_ligand++,
    smiles: '',
    atoms: [],
    bonds: [],
  }
  // 判断存在形式
  if (
    tmpSmile.smiles !== ''
    && (!JSON.stringify(mol4enum.ligand).includes(JSON.stringify(tmpSmile.smiles))
      || !JSON.stringify(mol4enum.ligand).includes(`"bonds":${JSON.stringify(tmpSmile.bonds)}`)
      || !JSON.stringify(mol4enum.ligand).includes(`"atoms":${JSON.stringify(tmpSmile.atoms)}`))
  ) {
    mol.smiles = tmpSmile.smiles
    mol.atoms = tmpSmile.atoms
    mol.bonds = tmpSmile.bonds
    mol4enum.ligand.push(JSON.parse(JSON.stringify(mol)))
    // tmpSmile.smiles=''
    // tmpSmile.atoms=[]
    // tmpSmile.bonds=[]
  }
  console.log(mol4enum.ligand)
}

function itemEdit(mol: any, index: number, molObjList: any[]) {
  molObjList[index].smiles = mol.smiles
  molObjList[index].atoms = mol.atoms
  molObjList[index].bonds = mol.bonds
  console.log('enum accept edit', mol)
}

watch(
  () => initMol.smiles,
  (val) => {
    tmpSmile.smiles = val
  },
)
</script>

<template>
  <div>
    <NPageHeader subtitle="分子结构枚举，子结构筛选 ( beta 1.0 )">
      <NGrid v-if="showAbout" :cols="5">
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
      <NDivider />
      <template #title>
        <NSpace>
          <a style="text-decoration: none; color: inherit">分子生成器</a>
          <NSwitch v-model:value="showLigandCore" size="medium" />
        </NSpace>
      </template>
      <template #avatar>
        <NAvatar
          color=" white"
          src="/favicon.svg"
          @click="showAbout = !showAbout"
        />
      </template>
      <template #extra>
        <NSpace>
          <NDropdown :options="options" placement="bottom-start">
            <NButton text>
              · · ·
            </NButton>
          </NDropdown>
        </NSpace>
      </template>
    </NPageHeader>
    <NCollapse :default-expanded-names="['1']" style="width: 60%">
      <NCollapseItem
        title=" RDKit位点标注"
        display-directive="show"
        class="collapse"
        name="1"
        style="width: 100%; height: 100%"
      >
        <NCard style="width: 100%; height: 100%">
          <NThing style="width: 100%; height: 100%">
            <template #description>
              <NInputGroup class="inputG">
                <NInput
                  v-model:value="inputText"
                  style="font-size: 20px; max-width: 90%"
                  type="text"
                  size="large"
                  clearable
                  placeholder="type smiles here"
                >
                  <template #suffix>
                    <NButton
                      color="#E29587"
                      circle
                      size="medium"
                      quaternary
                      @click="showModal = true"
                    >
                      <div class="i-carbon-cloud-satellite" />
                    </NButton>
                    <!-- modal画板区域 -->
                    <NModal v-model:show="showModal" display-directive="show">
                      <NCard
                        style="width: 1300px; height: 96vh"
                        :bordered="false"
                        size="huge"
                        role="dialog"
                        aria-modal="true"
                      >
                        <template #default>
                          <init-ketcher
                            @update-smiles="(smiles) => (inputText = smiles)"
                          />
                        </template>
                        <template #footer>
                          <NSpace justify="center">
                            <NButton @click="showModal = false">
                              取消
                            </NButton>
                            <NButton @click="showModal = false">
                              确定
                            </NButton>
                          </NSpace>
                        </template>
                      </NCard>
                    </NModal>
                  <!-- modal画板区域结束 -->
                  </template>
                </NInput>
                <NButton
                  size="large"
                  style="font-size: 20px; min-width: 10%"
                  color="#c471ed"
                  @click="drawMol"
                >
                  <template #icon>
                    <NIcon>
                      <div class="i-ion-color-palette-outline" />
                    </NIcon>
                  </template>
                  绘制
                </NButton>
              </NInputGroup>
            </template>
            <template #default>
              <div
                style="
                display: flex;
                align-items: center;
                justify-content: center;
              "
              >
                <!-- rdkit-sub v-bind="initMol"/ -->
                <editable-rdkit
                  :key="editRdkitKey"
                  v-bind="initMol"
                  style="width: 60%"
                  @update-mol="acceptMol"
                />
              </div>
            </template>
            <template #footer>
              <NSpace justify="space-between" style="width: 100%">
                <NButton
                  quaternary
                  :focusable="false"
                  @click="copy(initMol.smiles)"
                >
                  <template #icon>
                    <NIcon>
                      <div class="i-carbon-copy-file" />
                    </NIcon>
                  </template>
                  <NGradientText
                    word-break="normal"
                    style="max-width: 600px"
                    gradient="linear-gradient(90deg, red 0%, green 50%, blue 100%)"
                  >
                    {{ initMol.smiles }}
                  </NGradientText>
                </NButton>
                <NSpace justify="end">
                  <NButton
                    size="medium"
                    strong
                    secondary
                    round
                    type="info"
                    @click="addCore"
                  >
                    添加主核
                  </NButton>
                  <NButton
                    size="medium"
                    strong
                    secondary
                    round
                    type="info"
                    @click="addLigand"
                  >
                    添加配体
                  </NButton>
                  <NButton
                    szie="large"
                    round
                    color="#ff69b4"
                    @click="enumMol"
                  >
                    <template #icon>
                      <NIcon :size="20">
                        <div class="i-tabler-brand-apple-arcade" />
                      </NIcon>
                    </template>
                    开始枚举
                  </NButton>
                </NSpace>
              </NSpace>
            </template>
          </NThing>
        </NCard>
      </NCollapseItem>
    </NCollapse>
    <NDivider />
    <NGrid v-show="showLigandCore" x-gap="40" y-gap="20" :cols="2">
      <NGridItem>
        <NCollapse :default-expanded-names="['1']">
          <NCollapseItem title="主核结构" display-directive="show" name="1">
            <NGrid :cols="4" x-gap="8" y-gap="8">
              <NGridItem v-for="(item, index) in mol4enum.core" :key="item.id">
                <card-rdkit
                  v-bind="item"
                  @item-deleted="mol4enum.core.splice(index, 1)"
                  @it-edit-save="itemEdit($event, index, mol4enum.core)"
                />
              </NGridItem>
            </NGrid>
          </NCollapseItem>
        </NCollapse>
      </NGridItem>
      <NGridItem>
        <NCollapse :default-expanded-names="['1']">
          <NCollapseItem title="配体结构" display-directive="show" name="1">
            <NGrid :cols="4" x-gap="8" y-gap="8">
              <NGridItem
                v-for="(item, index) in mol4enum.ligand"
                :key="item.id"
              >
                <card-rdkit
                  v-bind="item"
                  @item-deleted="mol4enum.ligand.splice(index, 1)"
                  @it-edit-save="itemEdit($event, index, mol4enum.ligand)"
                />
              </NGridItem>
            </NGrid>
          </NCollapseItem>
        </NCollapse>
      </NGridItem>
    </NGrid>
    <NDivider />
    <NCollapse :default-expanded-names="['1']">
      <NCollapseItem display-directive="if" name="1">
        <template #header>
          <div style="background-color: #ddffff;padding: 14px;border-left: 6px solid #ccc;border-color: #2196F3; width:100%">
            结果输出
          </div>
        </template>
        <div style="padding-left:20px;padding-right:20px">
          <NScrollbar style="height:80vh;border: 1px solid #ccc;padding:5px;border-radius: 10px;width:100%;">
            <!-- <n-grid :cols="8" x-gap="8" y-gap="8">
        <n-grid-item v-for="(item, index) in resultData?.data.message" :key="index">
          <card-rdkit :smiles="item" style="width: 95%; height: 95%" />
        </n-grid-item>
      </n-grid> -->
            <Grid
              v-if="resultData?.data.message"
              :length="VLength"
              :page-size="VPageSize"
              :page-provider="pageProvider"
              :page-provider-debounce-time="1000"
              class="grid"
            >
              <template #probe>
                <div class="item" style="width:100%; padding-top:100%" />
              </template>

              <template #placeholder="{ style }">
                <div class="item" :style="style">
                  loanding...
                </div>
              </template>

              <template #default="{ item, style }">
                <div class="item" :style="style">
                  <card-rdkit :smiles="item" />
                </div>
              </template>
            </Grid>
          </NScrollbar>
        </div>
      </NCollapseItem>
    </NCollapse>
  </div>
</template>

<style scoped>
.collapse {
  padding: 20px;
  width: 100%;
}
.inputG {
  width: 100%;
  position: relative;
  border-radius: 3px 3px 3px 3px;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0px 1px 1px 0px rgba(0, 0, 0, 0.1);
}
.grid {
  display: grid;
  grid-gap: 8px;
  grid-template-columns: repeat(2, 1fr);
}

/* @media (min-width: 768px) {
  .grid {
    grid-template-columns: repeat(3, 1fr);
  }
} */
@media (min-width: 992px) {
  .grid {
    grid-template-columns: repeat(4, 1fr);
  }
}
@media (min-width: 1280px) {
  .grid {
    grid-template-columns: repeat(5, 1fr);
  }
}
@media (min-width: 1440px) {
  .grid {
    grid-template-columns: repeat(6, 1fr);
  }
}
@media (min-width: 1650px) {
  .grid {
    grid-template-columns: repeat(8, 1fr);
  }
}
@media (min-width: 2200px) {
  .grid {
    grid-template-columns: repeat(10, 1fr);
  }
}

.item {
  padding: 5px;
  text-align: center;
}
</style>
