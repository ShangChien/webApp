<script setup lang="ts">
import {
  NSpace,
  NButtonGroup,
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
         CopyFile, 
         Carbon,
         Move,
         Close,
         Minimize,
         ColorPalette
} from "@vicons/carbon";
import { useClipboard } from "@vueuse/core";
import initKetcher from "@/components/initKetcher.vue";
import classSites from "@/components/rdkitComponent/classSites.vue";
import tagMols from "@/components/rdkitComponent/tagMols.vue"
import editableRdkit from "@/components/rdkitComponent/editableRdkit.vue";
import { reactive, ref, inject, onMounted,computed } from "vue";
import type { Ref } from "vue";
import type { molData } from "@/components/types";
import { useDraggable,useElementSize } from '@vueuse/core'
import { useEditState } from '@/stores/editMol'
const EditStore = useEditState()

const initMol: molData = reactive({
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  qsmiles: "*~*",
  atoms: {1:[2,3]},
  bonds: {7:[3,4]},
  labels: [],
})
//获取浮动编辑框子组件的方法
const editMol=ref()
const isMounted = ref(false)
const undo = computed(() => isMounted ? editMol.value?.canUndo : false)
const redo = computed(() => isMounted ? editMol.value?.canRedo : false)
//刷新edit组件内部状态
const refreshKey = ref(1);
//定位浮动的位置
const el1 = ref<HTMLElement | null>(null)
const controlPin=ref<HTMLElement | null>(null)
const NBG=ref<HTMLElement | null>(null)
const { x, y, style } = useDraggable(el1, {
  initialValue: { x: 74, y: 9 },
})
const { width:widthBOX } = useElementSize(controlPin)
const { width:widthNBG } = useElementSize(NBG)
const mini=ref(true)
const visiualBox=inject('visiualBox')

const tmpSmile = editableRdkit.highlightMap
// const tmpSmile: molData = reactive({
//   smiles: initMol.smiles,
//   atoms: {},
//   bonds: {},
// });
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
onMounted(()=>{
  isMounted.value=true
})
</script>

<template>
<div style="width:100%;z-index:1;position: fixed" v-show="visiualBox">
  <div ref="el1" style="position: fixed;z-index:1;cursor:move" :style="style">
    <div v-if="!mini" :style="{'width':widthBOX-widthNBG-60+'px'}">
      <n-button
        size="tiny" color="#FFA48D" circle>
        <n-icon><move /></n-icon>
      </n-button>
    </div>
    <n-button v-else
              class="simplify" 
              style="position: fixed;cursor:grab"
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
       class="entireBox" 
       :style="{'left':x-60+'px','top':y-5+'px'}" 
       style="position:fixed"> 
    <div ref="controlPin" style="padding-bottom:5px;z-index:1">
      <n-button style="margin-left:0px;margin-right:5px;z-index:1" @click="mini=!mini" size="tiny" color="#7CBD99" circle>
        <n-icon><minimize /></n-icon>
      </n-button>
      <n-button style="z-index:1" @click="visiualBox=false" size="tiny" color="#F39BBA" circle>
        <n-icon ><close /></n-icon>
      </n-button>
    </div>
    <div style="position:absolute;top:2px;right:4px" ref="NBG">
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
            style="font-size: 20px;max-width: 100%"
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
          <n-button color="#c471ed" @click="drawMol" >
            <n-icon :size="30">
              <color-palette />
            </n-icon>
          </n-button>
        </n-input-group>
      </template>
      <template #default>
      <div style="
              display: flex;
              align-items: center;
              justify-content: center;
              margin-top: -8px;
            " >
        <editable-rdkit
          ref="editMol"
          :key="refreshKey"
          v-bind="initMol"
          style="width: 100%; border-radius: 5px"
        />
      </div>
      </template>
      <template #footer>
      <div :style="{'width':widthBOX+'px'}"
           style="margin-top:-8px">
        <div class="i-fluent-emoji-label text-3xl" style="float:left; width:10%" />
        <div style="float:left; width:4%" class="text-xl" >: </div>
        <div style="float:left; width:86%"><tag-mols /></div>
      </div>
      <n-space justify="space-between" style="width:100%; padding-top:2px" >
        <class-sites style="margin-top: -4px;" />
        <n-button
          size="small"
          strong
          round
          type="info"
          @click="EditStore.addMol(tmpSmile)"
          >添加结构</n-button>  
      </n-space>
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
  width: 50px; 
  height: auto; 
  aspect-ratio: 1;
  background-image: linear-gradient(to top, #f3e7e9 0%, #e3eeff 99%, #e3eeff 100%);
  border-radius: 10px;
  box-shadow: 2px 2px 5px 3px rgba(0, 0, 0, 0.1);
}
</style>