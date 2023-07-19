<script setup lang="ts">
import { NButton, NInput, NInputGroup } from 'naive-ui'
import { onBeforeUnmount, onMounted, ref } from 'vue'

const smiles = defineModel<string>('smiles')
const mounted = ref<boolean>(false)
const src = ref<string>()
const refketcher = ref<HTMLIFrameElement | any>()
const iframeWin = ref<any>(null)
function handleMessage(event: {
  data: { cmd: any; params: { data: null } }
}) {
  // console.log(event)
  setTimeout(() => mounted.value = true, 210)
  switch (event.data.cmd) {
    case 'postSmiles':
      smiles.value = event.data.params.data
      break
  }
}

function sendMessage() {
  iframeWin.value.postMessage(
    {
      cmd: 'setMole',
      params: {
        smiles: smiles.value,
      },
    },
    '*',
  )
}
function getMessage() {
  iframeWin.value.postMessage(
    {
      cmd: 'getSmiles',
      params: undefined,
    },
    '*',
  )
}
onMounted(() => {
  src.value = 'static/ketcher/index.html'
  window.addEventListener('message', handleMessage)
  iframeWin.value = refketcher.value.contentWindow
  // console.log(iframeWin.value)
})
onBeforeUnmount(() => {
  window.removeEventListener('message', handleMessage)
})
defineExpose({
  sendMessage,
  mounted,
})
</script>

<template>
  <div class="flex flex-col flex-nowrap justify-start items-center gap-1em">
    <NInputGroup class="flex-none">
      <NButton
        size="large"
        style="font-size: 20px"
        color="#6190E8"
        @click="getMessage()"
      >
        Get
      </NButton>
      <NInput
        v-model:value="smiles"
        style="font-size: 20px"
        class="ninput"
        type="text"
        size="large"
        color="#e0c3fc"
        clearable
        placeholder="type smiles here"
      />
      <NButton
        size="large"
        style="font-size: 20px"
        color="#FFAFBD"
        @click="sendMessage()"
      >
        Draw
      </NButton>
    </NInputGroup>
    <div class="flex-auto w-full flex justify-center items-center">
      <iframe
        v-show="mounted"
        id="ketcher"
        ref="refketcher"
        class="w-full h-full"
        frameborder="0"
        :src="src"
      />
      <div v-show="!mounted" class="i-svg-spinners-blocks-wave text-4xl bg-blue-3" />
    </div>
  </div>
</template>

<style>
</style>
