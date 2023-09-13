/// <reference types="vite/client" />
declare module "@/assets/RDKit_minimal.js"
declare module "pyodide/pyodide.js"
declare module 'naive-ui'
declare module "monaco-editor-vue3";
declare module 'pyodide/pyodide.mjs'
declare module "ngl/dist/ngl.js";
declare module 'svgo/lib/svgo.js';
declare module 'splitpanes';
declare interface Window {
  RDKit: any;
  py:any;
  MonacoEnvironment:any
}
declare module "*.vue" {
  import { defineComponent } from "vue";
  const Component: ReturnType<typeof defineComponent>;
  export default Component;
}
declare module 'virtual:pwa-register/vue' {
  // @ts-expect-error ignore when vue is not installed
  import type { Ref } from 'vue'

  export interface RegisterSWOptions {
    immediate?: boolean
    onNeedRefresh?: () => void
    onOfflineReady?: () => void
    onRegistered?: (registration: ServiceWorkerRegistration | undefined) => void
    onRegisterError?: (error: any) => void
  }

  export function useRegisterSW(options?: RegisterSWOptions): {
    needRefresh: Ref<boolean>
    offlineReady: Ref<boolean>
    updateServiceWorker: (reloadPage?: boolean) => Promise<void>
  }
}