import { createApp } from 'vue'
import { createPinia } from 'pinia'
import jquery from 'jquery'
Object.assign(window, { $: jquery})

import App from './App.vue'
import router from './router'

const app = createApp(App)

app.use(createPinia())
app.use(router)
app.mount('#app')
