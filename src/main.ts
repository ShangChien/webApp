import { createApp,ref,readonly } from "vue";
import { createPinia } from "pinia";
import App from "./App.vue";
import router from "./router";
import initRDKit from "@/components/rdkitComponent/RDKit";
//init SharedWorker(rdkit)
const myWorker = new SharedWorker(new URL('./worker/sharedWorker.js',import.meta.url))
myWorker.port.postMessage(null);
const app = createApp(App);

const rdkit = ref<any>();
initRDKit.then((res) => {
  rdkit.value = res; 
  app.provide('rdkit', readonly(rdkit.value))
  app.provide('myWorker', readonly(myWorker))
})



app.use(createPinia());
app.use(router);
app.mount("#app")
