<script setup lang="ts">
import { ref, h, reactive, onMounted, inject,computed } from "vue";
import type { Component, Ref } from "vue";
import type { molData,pgData } from "@/components/types";
import svgRdkit from "@/components/rdkitComponent/svgRdkit.vue";
import { useCopy } from '@/components/rdkitComponent/composable/useCopy'
import {
  NCard,
  NCheckbox,
  NTag,
  NDropdown,
  NButton,
  NIcon,
  NPopover,
} from "naive-ui";
import { useClipboard } from "@vueuse/core";
import { Dots } from "@vicons/tabler";
import { Edit, Delete, CopyFile } from "@vicons/carbon";
import { useEnumStore } from '@/stores/enumStore'

const emit = defineEmits(["itemChecked"]);
const props = defineProps<pgData>();
const currentEdit:Ref<{id:number;state:number}>=inject('currentEdit')
const enumStore = useEnumStore()
const editState=computed(()=>{
  if (currentEdit.value.state==2 && currentEdit.value.id==props.id) {
    return 'editing'
  } else if (currentEdit.value.state==1 && currentEdit.value.id==props.id) {
    return 'preview'
  } else {
    return 'normal'
  }
})

//可视加载组件
const checked = ref(false);
const cardView=ref()
const target = ref(null)

const { copy } = useClipboard();
const copytext = ref<string>(props.smiles);

const renderIcon = (icon: Component) => {
  return () => {
    return h(NIcon, null, {
      default: () => h(icon),
    });
  };
};
const options = [
  { label: "edit", key: "edit", icon: renderIcon(Edit) },
  { label: "copy", key: "copy", icon: renderIcon(CopyFile) },
  { label: "delete", key: "delete", icon: renderIcon(Delete) },
];
function handleDropOption(key: string) {
  console.log(key)
  switch (key) {
    case 'edit':
      currentEdit.value.id = props.id
      currentEdit.value.state = 1
      //console.log('editing',editState.value)
      break
    case 'copy':
      useCopy(``)
      //copy(copytext.value);
      console.log('copied')
      break
    case 'delete':
      enumStore.rmMolById(props.id)
      //console.log('delete ',props.id)
      break
    default:
      console.log('no option matched ')
  }
}
function quickView(){
  if (currentEdit.value.id != props.id) {
    currentEdit.value.id = props.id
    currentEdit.value.state = 1
    console.log('preview',editState.value)
  } else {
    currentEdit.value.id = 0
    currentEdit.value.state = 0
    console.log('exit preview')
  }
}
onMounted(() => {
  //console.log(props);
});
</script>
<template>
<div ref="target" class="w-100% h-100% ">
  <n-card ref="cardView" hoverable  >
    <template  #cover>
      <div class="ml-1 z-1 w-100%  flex justify-between">
        <n-checkbox v-model:checked="checked"  />
        <div>
          <n-tag  v-if="editState=='editing'" :bordered="false" type="warning" size="tiny">
            {{'editing:'+props.id}}
          </n-tag>
          <n-tag v-else-if="editState=='preview'" :bordered="false" type="primary" size="tiny">
            {{'preview:'+props.id}}
          </n-tag>
          <n-tag v-else :bordered="false" type="info" size="tiny">
            {{'id: '+props.id}}
          </n-tag>
        </div>
        <n-dropdown
          trigger="click"
          :options="options"
          placement="bottom-end"
          @select="handleDropOption">
          <n-button tertiary circle size="tiny" class="scale-90 mr-1">
            <template #icon>
              <n-icon><Dots /></n-icon>
            </template>
          </n-button>
        </n-dropdown>
      </div>
      <n-popover
        trigger="hover"
        placement="right"
        display-directive="if"
        to=".n-scrollbar"
      >
        <template #trigger>
          <div class="w-100% h-100%" @dblclick="quickView">
            <svg-rdkit v-bind="props" />
          </div> 
        </template>
        <span @click="copy(props.name_calc)"
          class="cursor-pointer rd-1 hover:bg-blue-50"
          >计算:{{ props.name_calc }}</span><br/>
        <span @click="copy(props.name_mat)"
          class="cursor-pointer rd-1 hover:bg-blue-50"
          >材料:{{ props.name_mat }}</span>
      </n-popover>
    </template>
  </n-card>
</div>
</template>
<style></style>
