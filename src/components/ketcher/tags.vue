<script setup lang='ts'>
import { computed,ref } from 'vue'
import { NButton,NCheckbox } from 'naive-ui'
import { useElementBounding,createReusableTemplate } from '@vueuse/core'

const props = defineProps<{tags:string[]}>()
const selected = defineModel<string[]>('selected')
const angle=ref(0)
const el_label = ref(null)
const expend  = ref(false)
const { height:labels_h,right:labels_r } = useElementBounding(el_label)
const [DefineTag, ReuseTag] = createReusableTemplate<{label:string }>()
const showExpand = computed(()=>{
  return (el_label.value?.lastElementChild.getBoundingClientRect().right > labels_r.value) 
    || (labels_h.value > 32) ? true : false})

async function select_all_or_none() {
  let len = await Promise.resolve(props.tags.length)
  if (len !== 0){
    selected.value = (selected.value.length!==0) ? [] : [...props.tags]
  }else{
    console.log('data items length', props.tags.length)
  }
}

async function reverse() {
  //翻转icon
  angle.value = (angle.value == 0) ? 180 : 0
  let arr= await Promise.resolve(props.tags)
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

</script>
<template>
<Transition name="fade">
  <div class="flex flex-nowrap justify-begin overflow-clip flex-auto items-start">
    <define-tag v-slot="{label}">
      <n-button class="m-0.5" round strong secondary size='tiny'
        :type="selected.includes(label) ? 'primary':'tertiary'"
        @click="select(label)">
        {{label}}
      </n-button>
    </define-tag>
    <div class="rd-1.5 bg-blue-100 h-24px
      flex-none flex flex-nowrap justify-around items-center 
      hover:bg-blue-200 ">
      <n-checkbox class=" ml-1.5 mr-1 " size="small"
      :checked="selected.length!==0 && selected.length===props.tags.length"
      :indeterminate="selected.length!==0 && selected.length!==props.tags.length"
      @click="select_all_or_none" >
      </n-checkbox>
      <div class="i-carbon-contrast m-1 ml-2 bg-indigo-400 cursor-pointer text-xl transition-210 hover:(bg-indigo-600)"
      :style="{ 'transform': `rotateY(${angle}deg)` }"
      @click="reverse"></div>
    </div>
    <div class="overflow-clip  flex-auto transition-height-210"
    :style="{height:labels_h+'px'}">
      <div ref="el_label" class="flex"
      :class="[expend ? 'flex-wrap':'flex-nowrap']">
        <reuse-tag v-for="tag in props.tags" :label="tag"></reuse-tag>
      </div>
    </div>
    <div v-if="showExpand" class="rd-2 ml-2 p-0.5 mt-0.5 bg-blue-100 flex-none hover:bg-blue-200">
      <div @click="()=>{expend = !expend;}"
      :class="[expend ? 'i-ion-chevron-collapse' : 'i-ion-chevron-expand']"
      class=" text-l p-1 c-zinc-400 hover:(c-slate-600 cursor-pointer)" />
    </div>
  </div>
</Transition>
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