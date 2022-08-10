import Components from "unplugin-vue-components/vite";
import { VitePWA } from 'vite-plugin-pwa'
import { defineConfig } from "vite";
import vue from "@vitejs/plugin-vue";
import vueJsx from "@vitejs/plugin-vue-jsx";
import VueTypeImports from "vite-plugin-vue-type-imports";
import mkcert from "vite-plugin-mkcert";
import {
  AntDesignVueResolver,
  NaiveUiResolver,
} from "unplugin-vue-components/resolvers";
import { resolve } from "path";

Components({
  resolvers: [AntDesignVueResolver(), NaiveUiResolver()],
});

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [
    vue(),
    mkcert(),
    vueJsx(),
    VueTypeImports(),
    VitePWA({
      includeAssets: ['favicon.svg', 'favicon.ico', 'robots.txt', 'apple-touch-icon.png'],
      manifest: {
        name: 'SUNERA App',
        short_name: 'App Space',
        description: 'sunera space components',
        theme_color: '#ffffff',
        icons: [
          {
            src: 'android-chrome-192x192.png',
            sizes: '192x192',
            type: 'image/png',
          },
          {
            src: 'android-chrome-384x384.png',
            sizes: '384x384',
            type: 'image/png',
          },
          {
            src: 'android-chrome-384x384.png',
            sizes: '384x384',
            type: 'image/png',
            purpose: 'any maskable',
          }
        ],
      },
      workbox: {
        maximumFileSizeToCacheInBytes: 10000000,
        globPatterns: ['**/*{js,css,html,ico,png,svg,pdb,sdf,wasm}']
      }
    }),
  ],
  server: {
    host: "0.0.0.0",
    https: true,
    cors: true, // 默认启用并允许任何源
    port: 8080,
    //反向代理配置，注意rewrite写法，开始没看文档在这里踩了坑
    proxy: {
      "^/api": {
        target: "http://192.168.2.160:8000", //代理接口
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, ""),
      },
    },
  },
  resolve: {
    alias: {
      "@": resolve(__dirname,"./src"),
    },
  },
});
