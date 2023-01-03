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
    nextId: 19,
  }),
  getters: {
    getByLabel(state) {
      return (label:string)=> state.mols.filter(
        (mol) => {
          return mol.labels.includes(label)
        }
      )
    },
    getById(state) {
      return (id:number)=> state.mols.filter(
        (mol) => mol.id === id
      )
    },
    getByType(state) {
      return (type:'core'|'ligand'|'mole')=> state.mols.filter(
        (mol) => mol.type == type
      )
    },
    getByTypeAndLabel(state) {
      return (type:'core'|'ligand'|'mole',label:string)=> state.mols.filter(
        (mol) => mol.type === type && label in mol?.labels
      )
    },
    getAllLabels(state):string[]|null[] {
      let init = []
      state.mols.forEach((item:molData)=>{
        init = init.concat(item.labels)
      })
      return Array.from(new Set(init))
    }
  },
  actions: {
    async addMol(x:molData) {
      Promise.resolve(x).then((res)=>{
        res.id=this.nextId++
        this.mols.push(res)
      }).catch(e=>console.log('add mol error:',e)) 
    },
    async updateMol(x:molData) {
      Promise.resolve(x).then((res:molData)=>{
        const index = this.mols.findIndex((mol) => mol.id == res.id)
        if (index!=-1) {
          this.mols[index]=res
        }
      }).catch((e)=>console.log(e,':id do not exist in enumMol'))
    },
    async rmMolById(id:number) {
      Promise.resolve(id).then((res)=>{
        this.mols = this.mols.filter((mol) => mol.id != res)
      }).catch((e)=>console.log(e,'error'))
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
  persistedState:{
    includePaths:['mols','nextId']
  } 
});