import { defineStore } from 'pinia'
import { clear, del, get, set } from 'idb-keyval'
import type { record, recordFull } from '@/components/types'

export const useSearchStore = defineStore('search', {
  state: () => ({
    records: [] as record[],
    nextId: 0,
  }),
  getters: {
    getById(state) {
      return async (id: number) => {
        let response: { id: number; smiles: string }[]
        await get(id).then((res) => {
          response = res
          console.log(`get res ${id} success`)
        }).catch(err => console.log(`get res ${id} failed! ${err}`))
        const out = state.records.find(mol => mol.id === id)
        return { ...out, res: response }
      }
    },
  },
  actions: {
    async addRecord(x: recordFull) {
      this.nextId++
      this.records.push({
        id: this.nextId,
        timestamp: x.timestamp,
        conditions: x.conditions,
        length: x.res.length,
      })
      await set(this.nextId, x.res).then(() => console.log(`save res ${this.nextId} success`))
        .catch(err => console.log(`save res ${this.nextId} failed! ${err}`))
    },
    async rmRecordById(id: number) {
      await del(id).then(() => console.log(`del res ${id} success`))
        .catch(err => console.log(`del res ${id} failed: ${err}`))
      this.records = this.records.filter((record: record) => record.id !== id)
    },
    async clearAll() {
      await clear().then(() => console.log('clear all success'))
        .catch(err => console.log(`clear all failed: ${err}`))
      this.records = []
      this.nextId = 0
    },
  },
})
