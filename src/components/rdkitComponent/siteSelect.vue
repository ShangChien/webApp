<script setup lang="ts">
import {
  NSpace,
  NButton,
  NIcon,
  NThing,
  NInputGroup,
  NInput,
  NModal,
  NCard,
  NScrollbar,
} from "naive-ui";
import { 
         CloudSatellite, 
         Carbon,
         Move,
         Close,
         Minimize,
         ColorPalette
} from "@vicons/carbon";
import { useClipboard,useWindowSize } from "@vueuse/core";
import initKetcher from "@/components/initKetcher.vue";
import classSites from "@/components/rdkitComponent/classSites.vue";
import tagMols from "@/components/rdkitComponent/tagMols.vue"
import editableRdkit from "@/components/rdkitComponent/editableRdkit.vue";
import { reactive, ref, inject, onMounted,computed } from "vue";
import type { molData } from "@/components/types";
import { useDraggable,useElementSize } from '@vueuse/core'
import { useEnumStore } from '@/stores/enumStore'
const enumStore = useEnumStore()
const initMol: molData = reactive({
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  qsmiles: "*~*",
  atoms: {1:[2,3]},
  bonds: {7:[3,4]},
  labels: [],
})
const molAddType = inject('molAdd',{molType:'ligand',info:'添加分子'})
//获取浮动编辑框子组件的方法
const editMol=ref()
const isMounted = ref(false)
const undo = computed(() => isMounted ? editMol.value?.canUndo : false)
const redo = computed(() => isMounted ? editMol.value?.canRedo : false)
const tmpSmile = computed(() => isMounted ? editMol.value.highlightMap: {})
//向pinia中添加分子，配体，主核
const storeAddMol = ()=>{
  tmpSmile.value['type'] = 'ligand'
  enumStore.addMol(tmpSmile.value)
  console.log(tmpSmile.value)
} 
//刷新edit组件内部状态
const refreshKey = ref(1);
//定位浮动的位置
const el1 = ref<HTMLElement | null>(null)
const controlPin=ref<HTMLElement | null>(null)
const NBG=ref<HTMLElement | null>(null)
const { width:windowWidth } = useWindowSize()
const { x, y, style } = useDraggable(el1, {
  initialValue: { x: windowWidth.value-325, y: 65 },
})
const { width:widthBOX } = useElementSize(controlPin)
const { width:widthNBG } = useElementSize(NBG)
const mini=ref(true)
const visiualBox=inject('visiualBox')

const drawMol=()=>{
  initMol.smiles = inputText.value;
  initMol.atoms = {};
  initMol.bonds = {};
  //console.log(initMol)
}
// const onReceiveMol=(mol:Ref<molData>)=>{
//   tmpSmile.smiles = mol.value.smiles;
//   tmpSmile.atoms = mol.value.atoms;
//   tmpSmile.bonds = mol.value.bonds;
//   //console.log(tmpSmile);
// }
const inputText = ref("");
const showModal = ref(false);
const { copy } = useClipboard();
function copySmile(){
  console.log(undo.value.value,redo.value.value)
  copy(initMol.smiles)
}

//分子类型
const currentType=ref<string>('mole')
const types=['ligand','core','mole']
onMounted(()=>{
  isMounted.value=true
})

</script>

<template>
<div class='z-3 fixed' v-show="visiualBox">
  <div ref="el1" class='z-3 fixed cursor-move' :style="style">
    <div v-if="!mini" :style="{'width':widthBOX-widthNBG-60+'px'}">
      <n-button
        size="tiny" color="#FFA48D" circle>
        <n-icon><move /></n-icon>
      </n-button>
    </div>
    <n-button v-else
              class="z-3 fixed cursor-grab simplify"
              :style="{left:x-60+'px',top:y-5+'px'}"
              @dblclick="mini=!mini"
              circle 
              color="#ff69b4" 
      >
      <template #icon>
        <n-icon :size="30">
          <carbon />
        </n-icon>
      </template>
    </n-button>
  </div>
  <div v-show="!mini" 
       class="entireBox fixed" 
       :style="{'left':x-60+'px','top':y-5+'px'}" > 
    <div ref="controlPin" class="pb-5px z-3" >
      <n-button style="margin-right:5px;z-index:3" @click="mini=!mini" size="tiny" color="#7CBD99" circle>
        <n-icon><minimize /></n-icon>
      </n-button>
      <n-button style="z-index:3" @click="visiualBox=false" size="tiny" color="#F39BBA" circle>
        <n-icon ><close /></n-icon>
      </n-button>
    </div>
    <div class='absolute top-2px right-4px' ref="NBG">
      <n-button circle :disabled="!undo" tertiary size="small" :focusable="false" @click="refreshKey++">
        <div class="i-ci-refresh-02 text-xl"></div> 
      </n-button>
      <n-button circle :disabled="!undo" tertiary size="small" :focusable="false" @click="editMol.undoRender()">
        <div class="i-flat-color-icons-undo text-xl"></div> 
      </n-button>
      <n-button circle :disabled="!redo" tertiary size="small" :focusable="false" @click="editMol.redoRender()">
        <div class="i-flat-color-icons-redo text-xl"></div> 
      </n-button>
      <n-button circle size="small" tertiary :focusable="false" @click="copySmile()">
        <div class="i-icon-park-copy text-xl"></div> 
      </n-button>
    </div>
    <n-scrollbar :x-scrollable="true" >
    <n-thing >
      <template #description>
        <n-input-group >
          <n-input
            v-model:value="inputText"
            class="text-20px max-w-100%"
            type="text"
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
                    <n-card class="w-1300px h-96vh" 
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
          <n-button color="#c471ed" @click="drawMol" >
            <n-icon :size="30">
              <color-palette />
            </n-icon>
          </n-button>
        </n-input-group>
      </template>
      <template #default>
      <div class="flex place-content-center mt--1" >
        <editable-rdkit
          ref="editMol"
          :key="refreshKey"
          v-bind="initMol"
          class="w-100% rd-5px"
        />
      </div>
      </template>
      <template #footer>
        <div class="gdview place-content-center mt--3">
          <class-sites />
          <div class="flex ma space-x-1 rd-2 bg-blue-900/20 p-1 " >
            <div v-for="(v,i) in types"
                    class="rd-1 text-sm p-1 leading-3 text-white 
                           focus:(outline-none ring-2) 
                           hover:(bg-white/[0.12] text-blue-400)"
                    :key="v"
                    :class="['tab-button shadow', { active: currentType === v }]"
                    @click="currentType= v" >
              {{ v }}
            </div>
          </div>
          <n-button size="small" type="info" secondary class="m-auto mr-1" @click="storeAddMol" >
            <template #icon> <div class="i-ion-plus-round text-teal-500" /></template>
            {{currentType}}
          </n-button>
        </div>
        <div :style="{'width':widthBOX+'px'}" >
          <div class="i-fluent-emoji-label text-3xl float-left mr-1 inline-block" />
          <div class="text-xl float-left inline-block" >:</div>
          <div class="text-xl float-left ml-2" :style="{'width':widthBOX-50+'px'}"><tag-mols /></div>
        </div>
      </template>
    </n-thing>
    </n-scrollbar>
  </div> 
</div>
</template>

<style scoped>
.entireBox{
  height: auto; 
  overflow:auto;
  resize:horizontal;
  background-image: linear-gradient(to top, #f3e7e9 0%, #e3eeff 99%, #e3eeff 100%);
  padding:5px;
  border-radius: 10px;
  box-shadow: 2px 2px 5px 5px rgba(0, 0, 0, 0.1);
}
.simplify{
  width: 40px; 
  height: auto; 
  aspect-ratio: 1;
  background-image: linear-gradient(to top, #f3e7e9 0%, #e3eeff 99%, #e3eeff 100%);
  border-radius: 10px;
  box-shadow: 2px 2px 5px 3px rgba(0, 0, 0, 0.1);
}
.gdview {
  display: grid;
  grid-template-columns: 90px 1fr 90px;
}
.tab-button.active {
  background-color: rgba(255, 255, 255, 1);
  color: rgba(37, 99, 235, 1)
}
</style>