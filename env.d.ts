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
