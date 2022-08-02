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
         CopyFile, 
         Carbon,
         Move,
         Close,
         Minimize,
         ColorPalette
} from "@vicons/carbon";
import { useClipboard } from "@vueuse/core";
import initKetcher from "@/components/initKetcher.vue";
import editableRdkit from "@/components/rdkitComponent/editableRdkit.vue";
import { reactive, ref, inject } from "vue";
import type { molData } from "@/components/types";
import { useDraggable,useElementSize } from '@vueuse/core'

const el1 = ref<HTMLElement | null>(null)
const controlPin=ref<HTMLElement | null>(null)
const { x, y, style } = useDraggable(el1, {
  initialValue: { x: 74, y: 9 },
})
const { width } = useElementSize(controlPin)

const mini=ref(true)
const visiualBox=inject('visiualBox')

const initMol: molData = reactive({
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  qsmiles: "*~*",
  atoms: [],
  bonds: [],
})
const tmpSmile = reactive({
  smiles: initMol.smiles,
  atoms: [],
  bonds: [],
});
const drawMol=()=>{
  initMol.smiles = inputText.value;
  tmpSmile.atoms = [];
  tmpSmile.bonds = [];
  //console.log(initMol)
}
const onReceiveMol=(mol: molData)=>{
  tmpSmile.smiles = mol.smiles;
  tmpSmile.atoms = mol.atoms;
  tmpSmile.bonds = mol.bonds;
  console.log(tmpSmile);
}
const inputText = ref("");
const showModal = ref(false);
const editRdkitKey = ref(1);
const { copy } = useClipboard();

</script>

<template>
<div style="width:100%;z-index:1;position: fixed" v-show="visiualBox">
  <div ref="el1" style="position: fixed;z-index:1" :style="style">
    <div v-if="!mini" :style="{'width':width*0.8+'px'}">
      <n-button
        size="tiny" color="#FFA48D" circle>
        <n-icon><move /></n-icon>
      </n-button>
    </div>
    <n-button v-else
              class="simplify" 
              style="position: fixed;"
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
       :style="{left:x-60+'px',top:y-5+'px'}" 
       style="position:fixed;"> 
    <div ref="controlPin" style="padding-bottom:5px;z-index:1">
      <n-button style="margin-left:0px;margin-right:5px;z-index:1" @click="mini=!mini" size="tiny" color="#7CBD99" circle>
        <n-icon><minimize /></n-icon>
      </n-button>
      <n-button style="z-index:1" @click="visiualBox=false" size="tiny" color="#F39BBA" circle>
        <n-icon ><close /></n-icon>
      </n-button>
    </div>
    <n-scrollbar :x-scrollable="true" style="resize:both">
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
          <n-button
            color="#c471ed"
            @click="drawMol"
          >
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
            " >
        <editable-rdkit
          :key="editRdkitKey"
          v-bind="initMol"
          @update-mol="onReceiveMol"
          style="width: 100%; border-radius: 5px"
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
            copy smiles
        </n-button>
          <n-button
            size="medium"
            strong
            secondary
            round
            type="info"
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