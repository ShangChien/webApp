import { defineStore } from "pinia";
import type { molData } from "@/components/types"

export const useEnumStore = defineStore('enum',{
  state: () => ({
    mols: [
      {id:1,smiles:"C1=CC1",labels:["三元环","芳环结构"],type:'ligand'},
      {id:2,smiles:"C1CCC1",labels:["四元环","烷环结构"],type:'ligand'},
      {id:3,smiles:"C1CCCC1",labels:["五元环","烷环结构"],type:'ligand'},
      {id:4,smiles:"C1CCCCC1",labels:["六元环","烷环结构"],type:'ligand'},
      {id:5,smiles:"C1CCCCCC1",labels:["七元环","烷环结构"],type:'ligand'},
      {id:6,smiles:"C1CCCCCCC1",labels:["八元环","烷环结构"],type:'ligand'},
      {id:7,smiles:"C1(C=CC=C2)=C2CC=C1",labels:["茚","芳环结构"],type:'ligand'},
      {id:8,smiles:"CC(C)C",labels:["叔丁基"],type:'ligand'},
      {id:9,smiles:"CC(C)(C)C1=CC(C(C)(C)C)=CC=C1",labels:["苯环两叔丁基","芳环结构"],type:'ligand'},
      {id:10,smiles:"C12=CC3C(C=CC=C3)=C1C=CC=C2",labels:["螺芴","芳环结构"],type:'ligand'},
      {id:11,smiles:"CC(C)C1=CC(C(C)C)=CC=C1",labels:["苯环二异丙基","芳环结构"],type:'ligand'},
      {id:12,smiles:"C12SC3C(C=CC=C3)=C1C=CC=C2",labels:["二苯并噻吩","芳环结构"],type:'ligand'},
      {id:13,smiles:"C12NC3C(C=CC=C3)=C1C=CC=C2",labels:["咔唑","芳环结构"],type:'ligand'},
      {id:14,smiles:"C12OC3C(C=CC=C3)=C1C=CC=C2",labels:["二苯并呋喃","芳环结构"],type:'ligand'},
      {id:15,smiles:"C12[Si]C3C(C=CC=C3)=C1C=CC=C2",labels:["硅螺芴","芳环结构"],type:'ligand'},
      {id:16,smiles:"C12SC3C(C=CC=C3)=C1C=CC=C2",labels:["二苯并噻吩","芳环结构"],type:'core'},
      {id:17,smiles:"C12NC3C(C=CC=C3)=C1C=CC=C2",labels:["咔唑","芳环结构"],type:'core'},
      {id:18,smiles:"C12OC3C(C=CC=C3)=C1C=CC=C2",labels:["二苯并呋喃","芳环结构"],type:'core'},
    ] as molData[],
    labels: [] as string[],
    nextId: 19,
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
    getByType(state) {
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
    async addMol(x:molData) {
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