import { resolve } from 'node:path'
import { defineConfig } from 'vite'
import { VitePWA } from 'vite-plugin-pwa'
import vue from '@vitejs/plugin-vue'
import vueJsx from '@vitejs/plugin-vue-jsx'
import mkcert from 'vite-plugin-mkcert'
import Unocss from 'unocss/vite'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [
    vue({
      script: {
        defineModel: true,
      },
    }),
    mkcert(),
    vueJsx(),
    Unocss(),
    VitePWA({
      strategies: 'generateSW',
      injectRegister: 'auto',
      includeAssets: ['/64.png', 'favicon.ico', 'robots.txt', '/180.png', 'favicon.svg'],
      manifest: {
        name: 'emm',
        short_name: 'emm',
        description: 'emm app',
        theme_color: '#ffffff',
        icons: [
          {
            src: '/192.png',
            sizes: '192x192',
            type: 'image/png',
          },
          {
            src: '/512.png',
            sizes: '512x512',
            type: 'image/png',
          },
          {
            src: '/512.png',
            sizes: '512x512',
            type: 'image/png',
            purpose: 'any maskable',
          },
        ],
        start_url: '/index.html',
      },
      workbox: {
        maximumFileSizeToCacheInBytes: 100000000,
        globPatterns: ['**/*.{js,css,html,ico,png,svg,pdb,sdf,wasm}'],
        navigateFallbackDenylist: [/^\/api/],
      },
    }),
  ],
  server: {
    host: '0.0.0.0',
    https: true,
    cors: true, // 默认启用并允许任何源
    port: 8080,
    // 反向代理配置，注意rewrite写法，开始没看文档在这里踩了坑
    proxy: {
      '/api/': {
        // target: 'http://api.vastlab.io:5050',
        target: 'http://192.168.2.233:5050', // 局域网代理接口
        // target: 'http://127.0.0.1:8000', // dev 本地接口
        changeOrigin: true,
        rewrite: (path: string) => path.replace(/^\/api/, ''),
      },
      // "^/src/worker/pyodide/.*\.whl$": {
      //   target: "https://cdn.jsdelivr.net/pyodide/v0.22.1/full", //代理接口
      //   changeOrigin: true,
      //   rewrite: (path:string) => {
      //     console.log(path.replace(/^\/src\/worker\/pyodide/, ''))
      //     return path.replace(/^\/src\/worker\/pyodide/, '')
      //   }
      // },
    },
  },
  resolve: {
    alias: {
      '@': resolve(__dirname, './src'),
    },
  },
})
