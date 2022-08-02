import { defineStore } from "pinia";
import type { mol4E } from "@/components/types"

export const useEnumCores = defineStore('enum',{
  state: () => ({
    cores: [] as mol4E[],
    nextId: 0,
  }),
  getters: {
    getCoreByLabel(state) {
      return (label:string)=> state.cores.filter(
        (mol) => mol.label ? (label in mol.label) : false
      )
    },
  },
  actions: {
    addMol(x:mol4E) {
      this.cores.push({
        id:this.nextId++,
        smiles:x.smiles,
        atoms:x.atoms,
        bonds:x.bonds,
        label:x.label,
      })
    },
    //delete mol by id
    rmMolById(id:number) {
      this.cores = this.cores.filter(
        (mol) => mol.id !== id
      )
    },
    //delete mol by label
    rmMolByLabel(label:string) {
      this.cores = this.cores.filter(
        (mol) => !(label in mol.label)
      )
    },
    //remove label in label array
    rmLabel(label:string) {
      for (var i = 0; i < this.cores.length; i++) {
        this.cores[i].label=this.cores[i].label.filter(
          (l) => l !== label
        ) 
      }
    }, 
  },
});