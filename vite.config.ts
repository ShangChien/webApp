import { fileURLToPath, URL } from "url";
import Components from "unplugin-vue-components/vite";
import { VitePWA } from 'vite-plugin-pwa'
import { defineConfig } from "vite";
import vue from "@vitejs/plugin-vue";
import vueJsx from "@vitejs/plugin-vue-jsx";
import VueTypeImports from "vite-plugin-vue-type-imports";
import mkcert from "vite-plugin-mkcert";
import UnoCSS from "unocss/vite";
import UnocssIcons from "@unocss/preset-icons";
import {
  AntDesignVueResolver,
  NaiveUiResolver,
} from "unplugin-vue-components/resolvers";

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
    UnoCSS({
      // 但 `presets` 被指定时，默认的预设将会被禁用，
      // 因此你可以在你原有的 App 上使用纯 CSS 图标而不需要担心 CSS 冲突的问题。
      presets: [
        UnocssIcons({
          // 其他选项
          prefix: "i-",
          extraProperties: {
            display: "inline-block",
          },
        }),
        // presetUno() - 取消注释以启用默认的预设
      ],
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
      "@": fileURLToPath(new URL("./src", import.meta.url)),
    },
  },
  base: "./",
});
