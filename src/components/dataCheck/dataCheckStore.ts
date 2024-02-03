import { ref } from 'vue'
import { defineStore } from 'pinia'

export const useDataCheckStore = defineStore('dataCheck', () => {
  const currentTab = ref<string>('files')
  const fileType = ref<string>('csv')
  const result = ref<any>('')

  return { currentTab, fileType, result }
})

export interface Spectrum {
  name: string
  nm: number
  intensity: number
}
export interface SpectrumFromDB {
  name: string
  raw_arr: number[][]
  peaks_arr: number[]
}
