<script setup lang='ts'>
import { computed,ref,h,watch,nextTick, onMounted } from 'vue'
import { NButton,NCheckbox,NPopover,NDynamicTags,NTag,NInput } from 'naive-ui'
import { useElementBounding,useElementSize,createReusableTemplate } from '@vueuse/core'

const props = defineProps<{tags:string[]}>()
const tags= ref(props.tags)
const selected = defineModel<string[]>('selected')
const angle=ref(0)
const el_label = ref(null)
const out_label = ref(null)
const showPopover = ref(false)
const clickOut=async()=>{
  setTimeout(()=>{showPopover.value=false},0)
}
const changeShowPopover=()=>showPopover.value=!showPopover.value

const expend  = ref(true)
const { right:labels_r } = useElementBounding(el_label)
const { height:labels_h } = useElementSize(el_label)
const [DefineTag, ReuseTag] = createReusableTemplate<{label:string }>()
const showExpand = computed(()=>{
  return (el_label.value?.lastElementChild.getBoundingClientRect().right > labels_r.value) 
    || (labels_h.value > 32) ? true : false})

async function select_all_or_none() {
  let len = await Promise.resolve(tags.value.length)
  if (len !== 0){
    selected.value = (selected.value.length!==0) ? [] : [...tags.value]
  }else{
    console.log('data items length', tags.value.length)
  }
}

async function reverse() {
  //翻转icon
  angle.value = (angle.value == 0) ? 180 : 0
  let arr= await Promise.resolve(tags.value)
  selected.value= arr.filter(tag => ! selected.value.includes(tag))
}

function select(label:string){
  if (selected.value.includes(label)) {
    selected.value.splice(selected.value.indexOf(label), 1)
    selected.value=[...selected.value]
  }else{
    selected.value.push(label)
    selected.value=[...selected.value]
  }
}
function delateTag(label:string){
  if (selected.value.includes(label)) {
    selected.value.splice(selected.value.indexOf(label), 1)
    selected.value=[...selected.value]
  }
}
//标签展示录入
const inputValue = ref('')
const pushTag=()=>{
  if (!tags.value.includes(inputValue.value)){
    tags.value.push(inputValue.value)
    inputValue.value=''
  }else{
    console.log(inputValue.value,'has already exist in tags')
  }
}
onMounted(()=>{

})

</script>
<template>
<n-popover trigger="manual" placement="bottom-start" width="trigger" raw class="rd-2 box-border bg-white"
  :show-arrow="false"
  :show="showPopover" 
  @clickoutside="clickOut">
  <template #trigger>
    <div class="flex flex-nowrap justify-around pr-1 pl-1 max-w-full w-full items-center">
      <div class="rd-2 m-0 p-0 h-22px w-28px bg-blue-200 cursor-pointer
        flex-none flex justify-center items-center hover:bg-blue-200"
        @click="changeShowPopover">
        <div class="i-material-symbols-settings-ethernet-rounded  m-0 p-0 c-blue-400 " />
      </div>
      <div class="overflow-hidden flex-auto transition-height-210 ml-1">
        <div ref="out_label" class="flex items-center flex-nowrap gap-1">
          <n-tag v-for="item in selected" type="success" round size="small" :closable="true"
          @close="delateTag(item)">
            {{item}}
          </n-tag>
          <span @click="changeShowPopover" v-if="selected.length==0" 
            class="w-full text-gray-4 text-center cursor-pointer">
            click to select tags
          </span>
        </div>
      </div>   
    </div>
  </template>
    <div class="flex flex-nowrap justify-begin max-w-full p-1 rd-2 box-border items-start relative">
      <define-tag v-slot="{label}">
        <n-button class="m-0.5" round strong secondary size='tiny'
          :type="selected.includes(label) ? 'primary':'tertiary'"
          @click="select(label)">
          {{label}}
        </n-button>
      </define-tag>
      <div class="flex-none flex flex-nowrap justify-around items-center 
        rd-1.5 bg-blue-100 h-24px hover:bg-blue-200 ">
        <n-checkbox class=" ml-1.5 mr-1 " size="small"
          :checked="selected.length!==0 && selected.length===tags.length"
          :indeterminate="selected.length!==0 && selected.length!==tags.length"
          @click="select_all_or_none" >
        </n-checkbox>
        <div class="i-carbon-contrast m-1 ml-2 bg-indigo-400 cursor-pointer text-xl transition-210 hover:(bg-indigo-600)"
        :style="{ 'transform': `rotateY(${angle}deg)` }"
        @click="reverse"></div>
      </div>
      <div class="flex-auto overflow-hidden flex flex-col justify-center items-start
        box-border relative transition-height-210">
        <div class="overflow-hidden box-border transition-height-210 max-w-full"
        :style="{height:labels_h+'px'}">
          <div ref="el_label" class="flex "
            :class="[expend ? 'flex-wrap':'flex-nowrap']">
            <reuse-tag v-for="tag in tags" :label="tag"></reuse-tag>
          </div>
        </div>
        <div class="max-w-full w-full box-border p-1">
          <n-input v-show="expend" clearable class="w-full"
            v-model:value="inputValue"
            @keyup.enter="pushTag"
            round size="tiny" placeholder="新增tags">   
          </n-input>
        </div>
      </div>
      <div v-if="showExpand" class="rd-2 ml-2 m-0 p-0 h-24px bg-blue-100 
        flex-none flex justify-center items-center hover:bg-blue-200">
        <div @click="expend = !expend"
          :class="[expend ? 'i-ion-chevron-collapse' : 'i-ion-chevron-expand']"
          class=" text-xl m-0 p-0 c-zinc-400 cursor-pointer hover:c-slate-600" />
      </div>
    </div>
</n-popover>
</template>
<style scoped>
.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.21s cubic-bezier(0.4, 0, 0.2, 1);
}
.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}
</style>