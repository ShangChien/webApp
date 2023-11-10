import { ref } from 'vue'
import { defineStore } from 'pinia'
import type { dataResults } from '@/components/types'

export const usemlStore = defineStore('machineLearning', () => {
  const currentTab = ref<string>('files')
  const fileType = ref<string>('*.mol')
  const task = ref<string>('regression')
  const models = ref<string[]>(['HOMO'])
  const isInferencing = ref<boolean>(false)
  const result = ref<dataResults[]>([])

  return { currentTab, fileType, task, models, result, isInferencing }
})
