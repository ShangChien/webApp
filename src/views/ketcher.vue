<script setup lang="ts">
import { NButton, NInput, NInputGroup, NSpace } from "naive-ui";
import { ref, onMounted, onBeforeUnmount } from "vue";

const src = ref<string>("/static/index.html?");
const smiles = ref<string | any>("");
const refketcher = ref<HTMLIFrameElement | any>();
const iframeWin = ref<any>(null);
const showketcher = ref<boolean>(true);

onMounted(() => {
  window.addEventListener("message", handleMessage);
  iframeWin.value = refketcher.value.contentWindow;
});
const handleMessage = (event: {
  data: { cmd: any; params: { data: null } };
}) => {
  switch (event.data.cmd) {
    case "postSmiles":
      smiles.value = event.data.params.data;
      break;
  }
};

function sendMessage() {
  iframeWin.value?.postMessage(
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

onBeforeUnmount(() => {
  window.removeEventListener("message", handleMessage);
});
</script>

<template>
  <div>
    <n-space vertical>
      <n-input-group class="inputG">
        <n-button
          size="large"
          style="font-size: 20px"
          color="#FFAFBD"
          @click="sendMessage"
        >
          Draw
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
          color="#3fada8"
          @click="getMessage"
        >
          Get
        </n-button>
        <n-button
          type="info"
          strong
          secondary
          round
          size="large"
          style="font-size: 20px"
          @click="showketcher = !showketcher"
          >{{ showketcher === false ? "打开画板" : "隐藏画板" }}
        </n-button>
      </n-input-group>
      <div class="ketcher" v-show="showketcher">
        <iframe
          id="ifKetcher"
          frameborder="0"
          ref="refketcher"
          :src="src"
          width="1050"
          height="780"
        ></iframe>
      </div>
    </n-space>
  </div>
</template>

<style scoped>
.ketcher {
  width: 850px;
  position: relative;
  height: 635px;
  padding: 5px 5px;
  box-sizing: border-box;
  border: 1px #ced8e4 solid;
  border-right: 0;
  border-radius: 5px 5px 5px 5px;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 1px 2px 2px 1px rgba(0, 0, 0, 0.1);
}
.ketcher > iframe {
  position: absolute;
  transform: scale(0.8);
  left: -100px;
  top: -70px;
}
.inputG {
  width: 850px;
  position: relative;
  border-radius: 3px 20px 20px 3px;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0px 1px 1px 0px rgba(0, 0, 0, 0.1);
}
</style>
