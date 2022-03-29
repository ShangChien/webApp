/// <reference types="vite/client" />
declare module '@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js'
declare module '3dmol/build/3Dmol-nojquery-min.js'
declare module '*.vue' {
    import { defineComponent } from 'vue'
    const Component: ReturnType<typeof defineComponent>
    export default Component
}
  