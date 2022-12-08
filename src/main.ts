import { createApp,ref,readonly } from "vue";
import { createPinia } from "pinia";
//import piniaPluginPersistedstate from 'pinia-plugin-persistedstate'
import { createPersistedStatePlugin } from 'pinia-plugin-persistedstate-2'
import localforage from "localforage"
import 'uno.css'
import App from "./App.vue";
import router from "./router";
import initRDKit from "@/components/rdkitComponent/RDKit";
//init SharedWorker(rdkit)
const myWorker = new SharedWorker(new URL('./worker/sharedWorker.js',import.meta.url),{name: 'vastLabSharedWorker',type: "module"})
myWorker.port.postMessage(null);

const siteType=ref<number>(0)
const molTags=ref<string[]>(["芳胺","咔唑","配体"])
const app = createApp(App);
app.provide('molTags',molTags)
app.provide('siteType',siteType)
const rdkit = ref<any>();
initRDKit.then((res) => {
  rdkit.value = res; 
  rdkit.value.prefer_coordgen(true);
  app.provide('rdkit', readonly(rdkit.value))
  app.provide('myWorker', readonly(myWorker))
})

const pinia = createPinia()
pinia.use(
  createPersistedStatePlugin({
    storage: {
      getItem: async (key) => {
        return localforage.getItem(key)
      },
      setItem: async (key, value) => {
        return localforage.setItem(key, value)
      },
      removeItem: async (key) => {
        return localforage.removeItem(key)
      },
    },
  }),
)
//pinia.use(piniaPluginPersistedstate)
app.use(pinia);
app.use(router);
app.mount("#app")
