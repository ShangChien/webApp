import { defineStore } from "pinia";
import type { condition } from "@/components/types"
import { set,get,del,clear } from 'idb-keyval';
interface records{
  id:number;
  timestamp:number;
  length:number;
  conditions:condition[];
}
interface recordFull extends records {
  res:{id:number,smiles:string}[];
}

export const useSearchStore = defineStore('search',{
  state: () => ({
    records: [] as records[],
    nextId: 1,
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
      await set(x.id,x.res).then(() => console.log(`save res ${x.id} success`))
      .catch((err) => console.log(`save res ${x.id} failed! ${err}`));
      this.records.push({
        id          : this.nextId++,
        timestamp   : x.timestamp,
        conditions  : x.conditions,
        length      : x.res.length
      })
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