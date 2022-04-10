import { fileURLToPath, URL } from 'url'
import Components from 'unplugin-vue-components/vite'
import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueJsx from '@vitejs/plugin-vue-jsx'
import VueTypeImports from 'vite-plugin-vue-type-imports'

import {
  AntDesignVueResolver,
  NaiveUiResolver
  
} from 'unplugin-vue-components/resolvers'
Components({
  resolvers: [
    AntDesignVueResolver(),
    NaiveUiResolver(),
  ],
})


// https://vitejs.dev/config/
export default defineConfig({
  plugins: [vue(), vueJsx(), VueTypeImports()],
  server: {
            host: '0.0.0.0'
          },
  resolve: {
    alias: {
      '@': fileURLToPath(new URL('./src', import.meta.url))
    }
  },
  //build: {
  //  rollupOptions: {
  //      external: ['jquery'],
  //      output: {
  //        globals: {
  //          jquery: '$'
  //        }
  //      }
  //  }
  //}
})
