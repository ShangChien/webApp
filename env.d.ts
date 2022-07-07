/// <reference types="vite/client" />
declare module "@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js";
declare module "ngl/dist/ngl.js";
declare interface Window {
  RDKit: any;
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