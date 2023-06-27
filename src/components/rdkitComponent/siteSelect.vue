<script setup lang="ts">
import {
  NButton,
  NCard,
  NIcon,
  NInput,
  NInputGroup,
  NModal,
  NScrollbar,
  NSpace,
  NThing,
} from 'naive-ui'
import {
  Carbon,
  Close,
  CloudSatellite,
  ColorPalette,
  Minimize,
  Move,
} from '@vicons/carbon'
import { useClipboard, useDraggable, useElementSize, useWindowSize } from '@vueuse/core'
import { computed, inject, onMounted, reactive, ref, watch, watchEffect } from 'vue'
import type { Ref } from 'vue'
import initKetcher from '@/components/ketcher/initKetcher.vue'
import classSites from '@/components/rdkitComponent/classSites.vue'
import tagMols from '@/components/rdkitComponent/tagMols.vue'
import editableRdkit from '@/components/rdkitComponent/editableRdkit.vue'
import type { molData } from '@/components/types'

import { useEnumStore } from '@/stores/enumStore'

const currentEdit: Ref<{ id: number;state: number }> = inject('currentEdit')
const enumStore = useEnumStore()
const initMol: molData = reactive({ qsmiles: '*~*' })
const id = computed(() => currentEdit.value.id)
// 获取浮动编辑框子组件的方法
const editMol = ref()
const isMounted = ref(false)
const undo = computed(() => isMounted.value ? editMol.value?.canUndo : false)
const redo = computed(() => isMounted.value ? editMol.value?.canRedo : false)
const tmpSmile = computed(() => isMounted.value ? editMol.value.highlightMap : {})
// 向pinia中添加分子，配体，主核
const types = ['ligand', 'core', 'mole']
const types4Store = ref<'core' | 'ligand' | 'mole'>()
const labels4Store = ref<string[]>([])
const buttonAnimate = ref(true)
async function updateMol() {
  const mol4Store = tmpSmile.value
  mol4Store.type = types4Store.value
  mol4Store.qsmiles = ''
  mol4Store.labels = labels4Store.value
  if (id.value === 0) { // 新增结构
    enumStore.addMol(mol4Store)
    currentEdit.value.state = 3
  }
  else { // 更新结构
    mol4Store.id = id.value
    await enumStore.updateMol(mol4Store)
    // change to complited state for updating checked-list
    currentEdit.value.state = 3
  }
  buttonAnimate.value = false
  setTimeout(() => {
    buttonAnimate.value = true
  }, 800)
  console.log(tmpSmile.value)
}
watchEffect(() => {
  const molInEnum = ref(enumStore.getById(id.value)[0])
  initMol.smiles = molInEnum.value?.smiles ?? ''
  initMol.atoms = molInEnum.value?.atoms ?? {}
  initMol.bonds = molInEnum.value?.bonds ?? {}
  initMol.labels = molInEnum.value?.labels ?? []
  initMol.type = molInEnum.value?.type ?? 'mole'
  types4Store.value = initMol.type
  labels4Store.value = initMol.labels
  // refreshKey.value++
  // console.log(molInEnum.value)
},
)
watch(
  [undo, types4Store, labels4Store],
  ([nUndo, nType, nLabel]) => {
    if (nUndo || (nType !== initMol.type) || (nLabel.toString() !== initMol.labels.toString()))
      currentEdit.value.state = 2
      // console.log('state editing:',currentEdit.value.state )

    else
      currentEdit.value.state = 1
      // console.log('state view:',currentEdit.value.state )
  },
  {
    flush: 'post',
  },
)
// 刷新edit组件内部状态
const refreshKey = ref(1)
// 定位浮动的位置
const el1 = ref<HTMLElement | null>(null)
const controlPin = ref<HTMLElement | null>(null)
const NBG = ref<HTMLElement | null>(null)
const { width: windowWidth } = useWindowSize()
const { x, y, style } = useDraggable(el1, {
  initialValue: { x: windowWidth.value - 400, y: 230 },
})
const { width: widthBOX } = useElementSize(controlPin)
const { width: widthNBG } = useElementSize(NBG)
const mini = ref(false)
const visiualBox = inject('visiualBox')

const inputText = ref('')
const showModal = ref(false)
const { copy } = useClipboard()
function copySmile() {
  console.log(undo.value.value, redo.value.value)
  copy(initMol.smiles)
}

function drawMol() {
  initMol.smiles = inputText.value
  initMol.atoms = {}
  initMol.bonds = {}
  // console.log(initMol)
}

onMounted(() => {
  isMounted.value = true
})
</script>

<template>
  <div v-show="visiualBox" class="z-3 fixed">
    <div ref="el1" class="z-3 fixed cursor-move" :style="style">
      <div v-if="!mini" :style="{ width: `${widthBOX - widthNBG - 60}px` }">
        <NButton
          size="tiny" color="#FFA48D" circle
        >
          <NIcon><Move /></NIcon>
        </NButton>
        {{ id === 0 ? 'mode: add new' : `mode: update-${id}` }}
      </div>
      <NButton
        v-else
        class="z-3 fixed cursor-grab simplify"
        :style="{ left: `${x - 60}px`, top: `${y - 5}px` }"
        circle
        color="#ff69b4"
        @dblclick="mini = !mini"
      >
        <template #icon>
          <NIcon :size="30">
            <Carbon />
          </NIcon>
        </template>
      </NButton>
    </div>
    <div
      v-show="!mini"
      class="entireBox fixed min-w-333px"
      :style="{ left: `${x - 60}px`, top: `${y - 5}px` }"
    >
      <div ref="controlPin" class="pb-5px z-3">
        <NButton style="margin-right:5px;z-index:3" size="tiny" color="#7CBD99" circle @click="mini = !mini">
          <NIcon><Minimize /></NIcon>
        </NButton>
        <NButton style="z-index:3" size="tiny" color="#F39BBA" circle @click="visiualBox = false">
          <NIcon><Close /></NIcon>
        </NButton>
      </div>
      <div ref="NBG" class="absolute top-2px right-4px">
        <NButton circle :disabled="!undo" tertiary size="small" :focusable="false" @click="refreshKey++">
          <div class="i-ic-twotone-settings-backup-restore text-xl" />
        </NButton>
        <NButton circle :disabled="!undo" tertiary size="small" :focusable="false" @click="editMol.undoRender()">
          <div class="i-flat-color-icons-undo text-xl" />
        </NButton>
        <NButton circle :disabled="!redo" tertiary size="small" :focusable="false" @click="editMol.redoRender()">
          <div class="i-flat-color-icons-redo text-xl" />
        </NButton>
        <NButton circle size="small" tertiary :focusable="false" @click="copySmile()">
          <div class="i-icon-park-copy text-xl " />
        </NButton>
      </div>
      <NScrollbar :x-scrollable="true">
        <NThing>
          <template #description>
            <NInputGroup>
              <NInput
                v-model:value="inputText"
                class="text-20px max-w-100%"
                type="text"
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
                    <NIcon
                      color="#D66D75"
                      size="30"
                      :component="CloudSatellite"
                    />
                  </NButton>
                  <!-- modal画板区域 -->
                  <NModal v-model:show="showModal" display-directive="show">
                    <NCard
                      class="w-1300px h-96vh"
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
              <NButton color="#c471ed" @click="drawMol">
                <NIcon :size="30">
                  <ColorPalette />
                </NIcon>
              </NButton>
            </NInputGroup>
          </template>
          <template #default>
            <div class="flex place-content-center mt--1">
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
              <div class="flex ma space-x-1 rd-2 bg-blue-900/20 p-1 ">
                <div
                  v-for="(v, _i) in types"
                  :key="v"
                  class="rd-1 text-sm p-1 leading-3 text-white
                           focus:(outline-none ring-2)
                           hover:(bg-white/[0.12] text-blue-400) tab-button shadow"
                  :class="[{ active: types4Store === v }]"
                  @click="(types4Store as string) = v"
                >
                  {{ v }}
                </div>
              </div>
              <NButton size="small" type="info" secondary class="m-auto mr-1" @click="updateMol">
                <div v-show="buttonAnimate">
                  <div v-if="id === 0" class="i-material-symbols-add-box-outline text-teal-500 text-2xl mr-1 ml--1" />
                  <div v-else class="i-fluent-save-sync-20-regular text-teal-500 text-2xl mr-1 ml--1" />
                </div>
                <div v-show="!buttonAnimate" class="i-carbon-checkmark-outline text-green  text-2xl mr-1 ml--1" />
                {{ types4Store }}
              </NButton>
            </div>
            <div :style="{ width: `${widthBOX}px` }">
              <div class="i-fluent-emoji-label text-3xl float-left mr-1 inline-block" />
              <div class="text-xl float-left inline-block">
                :
              </div>
              <div class="text-xl float-left ml-2" :style="{ width: `${widthBOX - 50}px` }">
                <tag-mols v-model:tags="labels4Store" />
              </div>
            </div>
          </template>
        </NThing>
      </NScrollbar>
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
