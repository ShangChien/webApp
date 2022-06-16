<script setup lang="ts">
import { ref, h, onMounted, reactive } from "vue";
import type { Component } from "vue";
import type { molData } from "@/components/types";
import svgRdkit from "@/components/svgRdkit.vue";
import editableRdkit from "@/components/editableRdkit.vue";
import {
  NCard,
  NCheckbox,
  NSpace,
  NDropdown,
  NButton,
  NIcon,
  NPopover,
  NModal,
} from "naive-ui";
import { useClipboard } from "@vueuse/core";
import { Dots } from "@vicons/tabler";
import { Edit, Delete, CopyFile } from "@vicons/carbon";
const emit = defineEmits(["itemChecked", "itEditSave", "itemDeleted"]);
const props = defineProps<molData>();

//handle Modal
const showModal = ref(false);
function onNegativeClick() {
  showModal.value = false;
}
const showOptions = ref(false);

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
  console.log(initSmile);
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

const checked = ref(false);
function visible() {
  if (checked.value == true) {
    showOptions.value = true;
  } else {
    showOptions.value = !showOptions.value;
  }
}

onMounted(() => {
  console.log(props);
});
</script>
<template>
  <n-card
    hoverable
    style="position: relative"
    @mouseover="visible"
    @mouseout="visible"
  >
    <template #cover>
      <n-space
        justify="space-between"
        v-show="showOptions"
        style="
          padding-left: 3%;
          padding-top: 2%;
          padding-right: 1%;
          position: absolute;
          z-index: 1;
          width: 95%;
          opacity: 1;
        "
      >
        <n-checkbox v-model:checked="checked" />
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
          <n-modal v-model:show="showModal" :mask-closable="false">
            <n-card
              style="width: 900px; height: 800px"
              title="位点重新选取"
              :bordered="false"
              size="huge"
              role="dialog"
              aria-modal="true"
            >
              <editable-rdkit
                v-bind="props"
                qsmiles="*~*"
                @update-mol="acceptMol"
                style="width: 70%"
              />
              <template #footer>
                <n-space>
                  <n-button @Click="onNegativeClick">取消</n-button>
                  <n-button @Click="onPositiveClick">确定</n-button>
                </n-space>
              </template>
            </n-card>
          </n-modal>
        </div>
      </n-space>
      <n-popover
        trigger="hover"
        placement="right"
        width="trigger"
        display-directive="if"
        :to="false"
        style="max-width: 100%"
      >
        <template #trigger>
          <svg-rdkit v-bind="props" style="width: 100%; height: 100%" />
        </template>
        <span style="word-break: break-word">
          {{ props.smiles }}<br />
          atoms:{{ props.atoms }}<br />
          bonds:{{ props.bonds }}
        </span>
      </n-popover>
    </template>
  </n-card>
</template>
<style></style>
