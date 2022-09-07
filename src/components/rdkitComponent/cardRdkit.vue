<script setup lang="ts">
import { ref, h, reactive,onMounted } from "vue";
import type { Component } from "vue";
import type { molData } from "@/components/types";
import svgRdkit from "@/components/rdkitComponent/svgRdkit.vue";
import {
  NCard,
  NSpin,
  NEmpty,
  NCheckbox,
  NSpace,
  NDropdown,
  NButton,
  NIcon,
  NPopover,
  NModal,
} from "naive-ui";
import { useClipboard,useMouseInElement } from "@vueuse/core";
import { Dots } from "@vicons/tabler";
import { Edit, Delete, CopyFile } from "@vicons/carbon";
import { defineAsyncComponent } from 'vue'
//异步加载mol编辑组件

const emit = defineEmits(["itemDeleted","itemChecked", "itEditSave"]);
const props = defineProps<molData>();
//可视加载组件
const checked = ref(false);
const cardView=ref()
const { isOutside }=useMouseInElement(cardView)
const target = ref(null)

//handle Modal
const showModal = ref(false);
function onNegativeClick() {
  showModal.value = false;
}

function onPositiveClick() {
  emit("itEditSave", initSmile);
  showModal.value = false;
}
const initSmile = reactive({
  smiles: props.smiles,
  atoms: ref(props.atoms),
  bonds: ref(props.bonds),
});
function acceptMol(mol: molData) {
  initSmile.smiles = mol.smiles;
  initSmile.atoms = mol.atoms;
  initSmile.bonds = mol.bonds;
  console.log(initSmile,"cardRdkit.vue acceptMol");
}

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
  if (key == "edit") {
    showModal.value = true;
  } else if (key == "copy") {
    copy(copytext.value);
  } else if (key == "delete") {
    emit("itemDeleted");
    // eslint-disable-next-line no-dupe-else-if
  } else if (key == "edit") {
    //
  }
}
onMounted(() => {
  console.log(props);
});
</script>
<template>
<div ref="target" style="width:100%;height:100%" >
  <!-- v-if="targetIsVisible" -->
  <n-card 
    hoverable
    ref="cardView"
  >
    <template #cover>
      <n-space
        justify="space-between"
        v-if="!isOutside||checked"
        style="
          padding-left: 3%;
          padding-top: 2%;
          padding-right: 1%;
          position: absolute;
          z-index: 1;
          width: 95%;
        "
      >
        <n-checkbox v-model:checked="checked" style="opacity: 0.8;" />
        <div>
          <n-dropdown
            trigger="click"
            :options="options"
            placement="bottom-end"
            @select="handleDropOption"
          >
            <n-button tertiary circle size="tiny">
              <template #icon>
                <n-icon><Dots /></n-icon>
              </template>
            </n-button>
          </n-dropdown>
          
        </div>
      </n-space>
      <n-popover
        trigger="hover"
        placement="right"
        width="trigger"
        display-directive="if"
        to=".n-scrollbar"
      >
        <template #trigger>
          <div style="width: 100%; height: 100%">
            <svg-rdkit  v-bind="props" />
          </div> 
        </template>
        <span style="word-break: break-word">
          {{ props.smiles }}<br />
          atoms:{{ props.atoms }}<br />
          bonds:{{ props.bonds }}
        </span>
      </n-popover>
    </template>
  </n-card>
</div>
</template>
<style></style>
