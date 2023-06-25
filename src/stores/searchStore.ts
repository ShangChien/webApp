import { defineStore } from "pinia";
import type { recordFull,records } from "@/components/types"
import { set,get,del,clear } from 'idb-keyval';

export const useSearchStore = defineStore('search',{
  state: () => ({
    records: [] as records[],
    nextId: 0,
  }),
  getters: {
    getById(state) {
      return async(id:number)=> {
        let response:{id:number,smiles:string}[]
        await get(id).then((res) => response=res)
        .then(() => console.log(`get res ${id} success`))
        .catch((err) => console.log(`get res ${id} failed! ${err}`));
        const out=state.records.filter((mol) => mol.id === id)
        return { ...out,res:response }
      }
    }
  },
  actions: {
    async addRecord(x:recordFull) {
      this.records.push({
        id          : this.nextId++,
        timestamp   : x.timestamp,
        conditions  : x.conditions,
        length      : x.res.length
      })
      await set(this.nextId,x.res).then(() => console.log(`save res ${this.nextId} success`))
      .catch((err) => console.log(`save res ${this.nextId} failed! ${err}`));
    },
    async rmRecordById(id:number) {
      await del(id).then(() => console.log(`del res ${id} success`))
      .catch((err) => console.log(`del res ${id} failed: ${err}`));
      this.records = this.records.filter(
        (mol) => mol.id !== id
      )
    },
    async clearAll() {
      await clear().then(() => console.log(`clear all success`))
      .catch((err) => console.log(`clear all failed: ${err}`));
      this.records = []
    }
  },
});