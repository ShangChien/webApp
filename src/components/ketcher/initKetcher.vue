<script setup lang="ts">
import { NButton, NInput, NInputGroup } from "naive-ui";
import { ref, onMounted, onBeforeUnmount,computed } from "vue";

const smiles = defineModel<string>('smiles')
const src = ref<string>();
const refketcher = ref<HTMLIFrameElement | any>();
const iframeWin = ref<any>(null);
const handleMessage = (event: {
  data: { cmd: any; params: { data: null } };
}) => {
  console.log(event)
  switch (event.data.cmd) {
    case "postSmiles":
      smiles.value = event.data.params.data;
      console.log(smiles.value)
      break;
  }
};

function sendMessage() {
  iframeWin.value.postMessage(
    {
      cmd: "setMole",
      params: {
        smiles: smiles.value,
      },
    },
    "*"
  );
}
function getMessage() {
  iframeWin.value.postMessage(
    {
      cmd: "getSmiles",
      params: void 0,
    },
    "*"
  );
}
onMounted(() => {
  src.value="static/ketcher/index.html"
  window.addEventListener("message", handleMessage);
  iframeWin.value = refketcher.value.contentWindow;
  //console.log(iframeWin.value)
});
onBeforeUnmount(() => {
  window.removeEventListener("message", handleMessage);
});
</script>

<template>
<div class="flex flex-col flex-nowrap justify-start items-center gap-1em">
  <n-input-group class="flex-none">
      <n-button
        size="large"
        style="font-size: 20px"
        color="#6190E8"
        @click="getMessage()"
      >
        Get
      </n-button>
      <n-input
        v-model:value="smiles"
        style="font-size: 20px"
        class="ninput"
        type="text"
        size="large"
        color="#e0c3fc"
        clearable
        placeholder="type smiles here"
      />
      <n-button
        size="large"
        style="font-size: 20px"
        color="#FFAFBD"
        @click="sendMessage()"
      >
        Draw
      </n-button>
  </n-input-group>
  <div class="flex-auto w-full ">
    <iframe
      class="w-full h-full"
      frameborder="0"
      ref="refketcher"
      :src="src"></iframe>
  </div>
</div>
</template>
<style>
</style>
