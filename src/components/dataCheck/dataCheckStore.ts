import { ref } from 'vue'
import { defineStore } from 'pinia'

export const useDataCheckStore = defineStore('dataCheck', () => {
  const currentTab = ref<string>('files')
  const fileType = ref<string>('csv')
  const result = ref<any>('')

  return { currentTab, fileType, result }
})

export function data2Spectrum(data: SpectrumFromDB): Spectrum[] {
  return data.raw_arr.map((item, index) => ({
    name: data.name,
    nm: 250 + index,
    intensity: item,
    peak: data.peaks_arr[401 - index] !== 0 ? data.peaks_arr[index] : null,
  }))
}

export interface Spectrum {
  name: string
  nm: number
  intensity: number
  peak?: number
}
export interface SpectrumFromDB {
  name: string
  raw_arr: number[]
  peaks_arr?: number[]
}
