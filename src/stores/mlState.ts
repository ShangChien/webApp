import { ref } from 'vue'
import type { dataResults } from '@/components/types'

// 全局状态，创建在模块作用域下
const fileType = ref('*.mol')
const task = ref('regression')
const models = ref<string[]>(['HOMO'])

const isInferencing = ref(false)
const result = ref<dataResults[]>([])

export function useMlState() {
  // 局部状态，每个组件都会创建
  // const localCount = ref(1)
  return { fileType, task, models, result, isInferencing }
}
