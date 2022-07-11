import { createApp,ref } from "vue";
import { createPinia } from "pinia";
import "uno.css";
import App from "./App.vue";
import router from "./router";
import initRDKit from "@/components/RDKit";

async () => {
  await initRDKit.then((res) => {
    window.RDKit = res;
    //console.log(res)
  })
}
const app = createApp(App);

app.use(createPinia());
app.use(router);
app.mount("#app")
