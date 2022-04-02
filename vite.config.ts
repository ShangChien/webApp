import { fileURLToPath, URL } from 'url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueJsx from '@vitejs/plugin-vue-jsx'
import VueTypeImports from 'vite-plugin-vue-type-imports'

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
