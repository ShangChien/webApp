<script setup lang="ts">
import { ref, h, reactive, onMounted, inject,computed } from "vue";
import type { Component, Ref } from "vue";
import type { molData } from "@/components/types";
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
const props = defineProps<molData>();
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
      useCopy(`<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 18.7,48.6 L 44.4,63.4 L 44.4,112.8 L 18.7,127.6 Z' style='fill:#9FACE6;fill-rule:evenodd;fill-opacity:1;stroke:#9FACE6;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-5 atom-0 atom-5' d='M 44.4,112.8 L 87.1,137.5 L 87.1,167.2 L 18.7,127.6 Z' style='fill:#9FACE6;fill-rule:evenodd;fill-opacity:1;stroke:#9FACE6;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-1 atom-1 atom-2' d='M 87.2,9.1 L 87.2,38.8 L 44.4,63.4 L 18.7,48.6 Z' style='fill:#9FACE6;fill-rule:evenodd;fill-opacity:1;stroke:#9FACE6;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-2 atom-2 atom-3' d='M 87.2,9.1 L 155.7,48.7 L 130.0,63.5 L 87.2,38.8 Z' style='fill:#9FACE6;fill-rule:evenodd;fill-opacity:1;stroke:#9FACE6;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-3 atom-3 atom-4' d='M 155.7,48.7 L 155.6,127.7 L 129.9,112.9 L 130.0,63.5 Z' style='fill:#9FACE6;fill-rule:evenodd;fill-opacity:1;stroke:#9FACE6;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-4 atom-4 atom-5' d='M 129.9,112.9 L 155.6,127.7 L 87.1,167.2 L 87.1,137.5 Z' style='fill:#9FACE6;fill-rule:evenodd;fill-opacity:1;stroke:#9FACE6;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<ellipse cx='31.5' cy='120.2' rx='19.3' ry='19.3' class='atom-0'  style='fill:#9FACE6;fill-rule:evenodd;stroke:#9FACE6;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<ellipse cx='31.6' cy='56.0' rx='19.3' ry='19.3' class='atom-1'  style='fill:#9FACE6;fill-rule:evenodd;stroke:#9FACE6;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<ellipse cx='87.2' cy='23.9' rx='19.3' ry='19.3' class='atom-2'  style='fill:#9FACE6;fill-rule:evenodd;stroke:#9FACE6;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<ellipse cx='142.8' cy='56.1' rx='19.3' ry='19.3' class='atom-3'  style='fill:#9FACE6;fill-rule:evenodd;stroke:#9FACE6;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<ellipse cx='142.8' cy='120.3' rx='19.3' ry='19.3' class='atom-4'  style='fill:#9FACE6;fill-rule:evenodd;stroke:#9FACE6;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<ellipse cx='87.1' cy='152.4' rx='19.3' ry='19.3' class='atom-5'  style='fill:#9FACE6;fill-rule:evenodd;stroke:#9FACE6;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 31.5,120.2 L 31.6,56.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 31.6,56.0 L 87.2,23.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 87.2,23.9 L 142.8,56.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 142.8,56.1 L 142.8,120.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 142.8,120.3 L 87.1,152.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 87.1,152.4 L 31.5,120.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 31.5,117.0 L 31.5,120.2 L 34.3,121.8' style='fill:none;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 31.6,59.2 L 31.6,56.0 L 34.3,54.4' style='fill:none;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 84.4,25.5 L 87.2,23.9 L 90.0,25.5' style='fill:none;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 140.0,54.5 L 142.8,56.1 L 142.8,59.3' style='fill:none;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 142.8,117.1 L 142.8,120.3 L 140.0,121.9' style='fill:none;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 89.9,150.8 L 87.1,152.4 L 84.3,150.8' style='fill:none;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
</svg>`)
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
    currentEdit.value.id =0
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
        width="trigger"
        display-directive="if"
        to=".n-scrollbar"
      >
        <template #trigger>
          <div class="w-100% h-100%" @dblclick="quickView">
            <svg-rdkit v-bind="props" />
          </div> 
        </template>
        <span style="word-break: break-word">
          id:{{ props.id }}<br />
          atoms:{{ props.atoms }}<br />
          bonds:{{ props.bonds }}
        </span>
      </n-popover>
    </template>
  </n-card>
</div>
</template>
<style></style>
