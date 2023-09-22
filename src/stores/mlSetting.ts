import { ref } from 'vue'

// 全局状态，创建在模块作用域下
const fileType = ref('*.mol(3D)')
const task = ref('regression')
const models = ref<string[]>(['HOMO'])

export function useMlSetting() {
  // 局部状态，每个组件都会创建
  // const localCount = ref(1)
  return { fileType, task, models }
}
