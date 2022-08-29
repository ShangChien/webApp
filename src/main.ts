import { createApp,ref,readonly } from "vue";
import { createPinia } from "pinia";
import 'uno.css'
import App from "./App.vue";
import router from "./router";
import initRDKit from "@/components/rdkitComponent/RDKit";
//init SharedWorker(rdkit)
const myWorker = new SharedWorker(new URL('./worker/sharedWorker.js',import.meta.url))
myWorker.port.postMessage(null);
const siteType=ref<number>(0)
const molTags=ref<string[]>(["label1","lab2","lab3","lab4"])
const app = createApp(App);
app.provide('molTags',molTags)
app.provide('siteType',siteType)
const rdkit = ref<any>();
initRDKit.then((res) => {
  rdkit.value = res; 
  app.provide('rdkit', readonly(rdkit.value))
  app.provide('myWorker', readonly(myWorker))
})



app.use(createPinia());
app.use(router);
app.mount("#app")
