import { defineStore } from "pinia";
import type { mol4E } from "@/components/types"

export const useEnumCores = defineStore('enum',{
  state: () => ({
    cores: [] as mol4E[],
    nextId: 0,
  }),
  getters: {
    getByLabel(state) {
      return (label:string)=> state.cores.filter(
        (mol) => label in mol.label
      )
    },
    getById(state) {
      return (id:number)=> state.cores.filter(
        (mol) => mol.id === id
      )
    }
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
    addLabelById(id:number,label:string) {
      const item = this.cores.find(
        (mol) => mol.id === id
      )
      if (item?.label.indexOf(label) === -1) {
        item.label.push(label)
      } else {
        console.log('label already exists') 
      }  
    },
    rmMolById(id:number) {
      this.cores = this.cores.filter(
        (mol) => mol.id !== id
      )
    },
    rmMolByLabel(label:string) {
      this.cores = this.cores.filter(
        (mol) => !(label in mol.label)
      )
    },
    rmLabel(label:string) {
      for (var i = 0; i < this.cores.length; i++) {
        this.cores[i].label=this.cores[i].label.filter(
          (le) => le !== label
        )
      }
    }, 
  },
});