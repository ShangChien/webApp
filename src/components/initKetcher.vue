<script setup lang='ts'>
import { NButton,NInput,NInputGroup,NSpace } from 'naive-ui'
import { ref,onMounted, onBeforeUnmount} from 'vue'

const src = ref<string>('/static/index.html?')
const smiles = ref<string|any>('')
const refketcher = ref<HTMLIFrameElement|any>()
const iframeWin = ref<any>(null)
const emit = defineEmits(['updateSmiles'])
onMounted(()=>{
  window.addEventListener("message", handleMessage)
  iframeWin.value=refketcher.value.contentWindow
  //console.log(iframeWin.value)
})
const handleMessage = (event: { data: { cmd: any; params: { data: null; }; }; }) => {
  switch (event.data.cmd) {
    case 'postSmiles':
      smiles.value= event.data.params.data
      emit('updateSmiles',smiles.value)
      break
  }
};

function sendMessage() {
  iframeWin.value?.postMessage({
    cmd: 'setMole',
    params: {
      smiles:smiles.value
    }
  },'*')
}
function getMessage() {
  iframeWin.value.postMessage({
    cmd: 'getSmiles',
    params: void 0
    },'*')
}


onBeforeUnmount(() => {
	window.removeEventListener("message", handleMessage);
});
</script>

<template>
  <n-space vertical style="position: relative;" >
    <n-input-group class="inputG" >
      <n-button size="large" 
                style="font-size:20px"
                color="#6190E8"
                @click="getMessage">
        Get
      </n-button>
      <n-input v-model:value="smiles"
               style="font-size:20px"
               class="ninput"
               type="text"
               size="large"
               color = "#e0c3fc"
               clearable
               placeholder="type smiles here" 
               />
      <n-button size="large"
                style="font-size:20px"
                color="#FFAFBD"
                @click="sendMessage">
         Draw
      </n-button>
    </n-input-group>
    <div class="ketcher">
      <iframe id="ifKetcher"
            frameborder="0"
            ref="refketcher"
            :src="src"
            
            ></iframe>
    </div>
  </n-space>
</template>

<style scoped>
.ketcher{
    width: 100%;
    height: 80vh;
    position: relative;
    padding: 5px 5px;
    box-sizing: border-box; 
    border: 1px #ced8e4 solid; 
    border-right: 0;
    border-radius: 5px 5px 5px 5px ; 
    box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 1px 2px 2px 1px  rgba(0, 0, 0, 0.1);
}
.ketcher>iframe {
  position: absolute;
  width:103.8%;
  height:83vh;
  left: -24px;
  top:-13px;
  transform: scale(0.96);
}
.inputG{
    width: 100%;
    position: relative;
    border-radius: 3px 3px 3px 3px;
    box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0px 1px 1px 0px rgba(0, 0, 0, 0.1);
}
</style>

