import { ref } from 'vue'
import { defineStore } from 'pinia'

export const useDataCheckStore = defineStore('dataCheck', () => {
  const currentTab = ref<string>('files')
  const fileType = ref<string>('csv')
  const result = ref<any>('')

  return { currentTab, fileType, result }
})

export interface spectrumData {
  name: string
  nm: number[]
  intensity: number[]
}
