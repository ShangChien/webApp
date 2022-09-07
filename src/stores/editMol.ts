import { defineStore } from "pinia";
import type { molData } from "@/components/types"

export const useEditState = defineStore('edit',{
  state: () => ({
    mols: [] as molData[],
    nextId: 0,
  }),
  getters: {
    //获取最后一个mol
    getLastMol(state) {
      return state.mols[state.mols.length-1]
    },
    getByLabel(state) {
      return (label:string)=> state.mols.filter(
        (mol) => label in mol?.[label]
      )
    },
    getById(state) {
      return (id:number)=> state.mols.filter(
        (mol) => mol.id === id
      )
    }
  },
  actions: {
    addMol(x:molData) {
      this.mols.push({
        id:this.nextId++,
        smiles:x.smiles,
        atoms:x.atoms,
        bonds:x.bonds,
        label:x.label,
      })
    },
    addLabelById(id:number,label:string) {
      const item = this.mols.find(
        (mol) => mol.id === id
      )
      if (item?.[label].indexOf(label) === -1) {
        item?.[label].push(label)
      } else {
        console.log('label already exists') 
      }  
    },
    rmMolById(id:number) {
      this.mols = this.mols.filter(
        (mol) => mol.id !== id
      )
    },
    rmMolByLabel(label:string) {
      this.mols = this.mols.filter(
        (mol) => !(label in mol?.[label])
      )
    },
    rmLabel(label:string) {
      for (var i = 0; i < this.mols.length; i++) {
        this.mols[i].label=this.mols[i]?.[label].filter(
          (le:string|null) => le !== label
        )
      }
    }, 
  },
});