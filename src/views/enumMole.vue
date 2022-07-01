<script setup lang="ts">
import {
  NSpace,
  NPageHeader,
  NGrid,
  NGridItem,
  NStatistic,
  NButton,
  NDropdown,
  NAvatar,
  NIcon,
  NSwitch,
  NCollapse,
  NThing,
  NCollapseItem,
  NInputGroup,
  NInput,
  NDivider,
  NModal,
  NCard,
  NGradientText,
} from "naive-ui";
import { ColorPaletteOutline } from "@vicons/ionicons5";
import { CloudSatellite, CopyFile } from "@vicons/carbon";
import { BrandAppleArcade } from "@vicons/tabler";
import { useClipboard } from "@vueuse/core";
import initKetcher from "../components/initKetcher.vue";
import editableRdkit from "@/components/editableRdkit.vue";
import cardRdkit from "@/components/cardRdkit.vue";
import { reactive, ref, watch } from "vue";
import type { molData } from "@/components/types";
import axios from "axios";
import { useAxios } from "@vueuse/integrations/useAxios";

const instance = axios.create({
  baseURL: "/api",
});
const { data, execute } = useAxios("/enum", 
                                   { method: "POST" }, 
                                   instance, 
                                   {immediate: false,}
                                   );
async function getData() {
  await execute(undefined, { data: mol4enum })
}

const inputText = ref("");
const showAbout = ref(false);
const showModal = ref(false);
const editRdkitKey = ref(1);
const showLigandCore = ref(true);
const options = [{ label: "Ketcher: 2.4.0" }, { label: "RDKit: 2022.3.2" }];
const initMol: molData = reactive({
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  qsmiles: "*~*",
  atoms: [],
  bonds: [],
});

function drawMol() {
  initMol.smiles = inputText.value;
  tmpSmile.atoms = [];
  tmpSmile.bonds = [];
  //console.log(initMol)
}

const { copy } = useClipboard();

const tmpSmile = reactive({
  smiles: initMol.smiles,
  atoms: [],
  bonds: [],
});
function acceptMol(mol: molData) {
  tmpSmile.smiles = mol.smiles;
  tmpSmile.atoms = mol.atoms;
  tmpSmile.bonds = mol.bonds;
  console.log(tmpSmile);
}

const mol4enum = reactive({
  core: [] as Array<molData>,
  ligand: [] as Array<molData>,
});
//addcore
let id_core = 1;
function addCore() {
  //editRdkitKey.value++ 是否刷新画板
  const mol: molData = {
    id: id_core++,
    smiles: "",
    atoms: [],
    bonds: [],
  };
  //判断存在形式
  if (
    tmpSmile.smiles != "" &&
    (JSON.stringify(mol4enum.core).indexOf(JSON.stringify(tmpSmile.smiles)) ==
      -1 ||
      JSON.stringify(mol4enum.core).indexOf(
        '"bonds":' + JSON.stringify(tmpSmile.bonds)
      ) == -1 ||
      JSON.stringify(mol4enum.core).indexOf(
        '"atoms":' + JSON.stringify(tmpSmile.atoms)
      ) == -1)
  ) {
    mol.smiles = tmpSmile.smiles;
    mol.atoms = tmpSmile.atoms;
    mol.bonds = tmpSmile.bonds;
    mol4enum.core.push(JSON.parse(JSON.stringify(mol)));
    //tmpSmile.smiles=''
    //tmpSmile.atoms=[]
    //tmpSmile.bonds=[]
  }
  console.log(mol4enum.core);
}
//addLigand
let id_ligand = 1;
function addLigand() {
  //editRdkitKey.value++
  const mol: molData = {
    id: id_ligand++,
    smiles: "",
    atoms: [],
    bonds: [],
  };
  //判断存在形式
  if (
    tmpSmile.smiles != "" &&
    (JSON.stringify(mol4enum.ligand).indexOf(JSON.stringify(tmpSmile.smiles)) ==
      -1 ||
      JSON.stringify(mol4enum.ligand).indexOf(
        '"bonds":' + JSON.stringify(tmpSmile.bonds)
      ) == -1 ||
      JSON.stringify(mol4enum.ligand).indexOf(
        '"atoms":' + JSON.stringify(tmpSmile.atoms)
      ) == -1)
  ) {
    mol.smiles = tmpSmile.smiles;
    mol.atoms = tmpSmile.atoms;
    mol.bonds = tmpSmile.bonds;
    mol4enum.ligand.push(JSON.parse(JSON.stringify(mol)));
    //tmpSmile.smiles=''
    //tmpSmile.atoms=[]
    //tmpSmile.bonds=[]
  }
  console.log(mol4enum.ligand);
}
// eslint-disable-next-line @typescript-eslint/no-explicit-any
function itemEdit(mol: any, index: number, molObjList: any[]) {
  molObjList[index].smiles = mol.smiles;
  molObjList[index].atoms = mol.atoms;
  molObjList[index].bonds = mol.bonds;
}

watch(
  () => initMol.smiles,
  (val) => {
    tmpSmile.smiles = val;
  }
);
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
        <a style="text-decoration: none; color: inherit">分子生成器</a>
        <n-switch v-model:value="showLigandCore" size="medium" />
      </n-space>
    </template>
    <template #avatar>
      <n-avatar
        @click="showAbout = !showAbout"
        color=" white"
        src="/favicon.svg"
      />
    </template>
    <template #extra>
      <n-space>
        <n-dropdown :options="options" placement="bottom-start">
          <n-button text>· · ·</n-button>
        </n-dropdown>
      </n-space>
    </template>
  </n-page-header>
  <n-collapse :default-expanded-names="['1']" style="width: 60%">
    <n-collapse-item
      title=" RDKit位点标注"
      display-directive="show"
      class="collapse"
      name="1"
      style="width: 100%; height: 100%"
    >
      <n-card style="width: 100%; height: 100%">
        <n-thing style="width: 100%; height: 100%">
          <template #description>
            <n-input-group class="inputG">
              <n-input
                v-model:value="inputText"
                style="font-size: 20px; max-width: 90%"
                type="text"
                size="large"
                clearable
                placeholder="type smiles here"
              >
                <template #suffix>
                  <n-button
                    @click="showModal = true"
                    color="#E29587"
                    circle
                    size="medium"
                    quaternary
                  >
                    <n-icon
                      color="#D66D75"
                      size="30"
                      :component="CloudSatellite"
                    />
                  </n-button>
                  <!--modal画板区域-->
                  <n-modal v-model:show="showModal" display-directive="show">
                    <n-card
                      style="width: 1300px; height: 96vh"
                      :bordered="false"
                      size="huge"
                      role="dialog"
                      aria-modal="true"
                    >
                      <template #default>
                        <init-ketcher
                          @update-smiles="(smiles) => (inputText = smiles)"
                        ></init-ketcher>
                      </template>
                      <template #footer>
                        <n-space justify="center">
                          <n-button @Click="showModal = false">取消</n-button>
                          <n-button @Click="showModal = false">确定</n-button>
                        </n-space>
                      </template>
                    </n-card>
                  </n-modal>
                  <!--modal画板区域结束-->
                </template>
              </n-input>
              <n-button
                size="large"
                style="font-size: 20px; min-width: 10%"
                color="#c471ed"
                @click="drawMol"
              >
                <template #icon>
                  <n-icon>
                    <ColorPaletteOutline />
                  </n-icon>
                </template>
                绘制
              </n-button>
            </n-input-group>
          </template>
          <template #default>
            <div
              style="
                display: flex;
                align-items: center;
                justify-content: center;
              "
            >
              <!--rdkit-sub v-bind="initMol"/-->
              <editable-rdkit
                :key="editRdkitKey"
                v-bind="initMol"
                @update-mol="acceptMol"
                style="width: 60%"
              />
            </div>
          </template>
          <template #footer>
            <n-space justify="space-between" style="width: 100%">
              <n-button
                quaternary
                :focusable="false"
                @click="copy(initMol.smiles)"
              >
                <template #icon>
                  <n-icon>
                    <CopyFile />
                  </n-icon>
                </template>
                <n-gradient-text
                  word-break="normal"
                  style="max-width: 600px"
                  gradient="linear-gradient(90deg, red 0%, green 50%, blue 100%)"
                >
                  {{ initMol.smiles }}
                </n-gradient-text>
              </n-button>
              <n-space justify="end">
                <n-button
                  size="medium"
                  strong
                  secondary
                  round
                  type="info"
                  @click="addCore"
                  >添加主核</n-button
                >
                <n-button
                  size="medium"
                  strong
                  secondary
                  round
                  type="info"
                  @click="addLigand"
                  >添加配体</n-button
                >
                <n-button
                  szie="large"
                  round
                  color="#ff69b4"
                  @click="getData"
                >
                  <template #icon>
                    <n-icon :size="20">
                      <BrandAppleArcade />
                    </n-icon>
                  </template>
                  开始枚举
                </n-button>
              </n-space>
            </n-space>
          </template>
        </n-thing>
      </n-card>
    </n-collapse-item>
  </n-collapse>
  <n-divider />
  <n-grid x-gap="40" y-gap="20" :cols="2" v-show="showLigandCore">
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="主核结构" display-directive="show" name="1">
          <n-grid :cols="4" x-gap="8" y-gap="8">
            <n-grid-item v-for="(item, index) in mol4enum.core" :key="item.id">
              <card-rdkit
                v-bind="item"
                @item-deleted="mol4enum.core.splice(index, 1)"
                @it-edit-save="itemEdit($event, index, mol4enum.core)"
              />
            </n-grid-item>
          </n-grid>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
    <n-grid-item>
      <n-collapse :default-expanded-names="['1']">
        <n-collapse-item title="配体结构" display-directive="show" name="1">
          <n-grid :cols="4" x-gap="8" y-gap="8">
            <n-grid-item
              v-for="(item, index) in mol4enum.ligand"
              :key="item.id"
            >
              <card-rdkit
                v-bind="item"
                @item-deleted="mol4enum.ligand.splice(index, 1)"
                @it-edit-save="itemEdit($event, index, mol4enum.ligand)"
              />
            </n-grid-item>
          </n-grid>
        </n-collapse-item>
      </n-collapse>
    </n-grid-item>
  </n-grid>
  <n-divider />
  <n-collapse :default-expanded-names="['1']">
    <n-collapse-item title="结果输出" display-directive="show" name="1">
      <n-grid :cols="8" x-gap="8" y-gap="8">
        <n-grid-item v-for="(item, index) in data?.message" :key="index">
          <card-rdkit :smiles="item" style="width: 95%; height: 95%" />
        </n-grid-item>
      </n-grid>
    </n-collapse-item>
  </n-collapse>
</template>

<style>
.collapse {
  padding: 5px;
  width: 100%;
}
.inputG {
  width: 100%;
  position: relative;
  border-radius: 3px 3px 3px 3px;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0px 1px 1px 0px rgba(0, 0, 0, 0.1);
}
</style>
