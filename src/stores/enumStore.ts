import { defineStore } from "pinia";
import type { molData } from "@/components/types"

export const useEnumStore = defineStore('enum',{
  state: () => ({
    mols: [] as molData[],
    nextId: 0,
  }),
  getters: {
    getByLabel(state) {
      return (label:string)=> state.mols.filter(
        (mol) => label in mol?.labels
      )
    },
    getById(state) {
      return (id:number)=> state.mols.filter(
        (mol) => mol.id === id
      )
    },
    getBytype(state) {
      return (type:'core'|'ligand'|'mole')=> state.mols.filter(
        (mol) => mol.type === type
      )
    },
    getByTypeAndLabel(state) {
      return (type:'core'|'ligand'|'mole',label:string)=> state.mols.filter(
        (mol) => mol.type === type && label in mol?.labels
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
        labels:x.labels,
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
        this.mols[i].labels=this.mols[i]?.labels.filter(
          (le:string|null) => le !== label
        )
      }
    }, 
  },
});