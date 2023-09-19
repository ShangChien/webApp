<script setup lang='ts'>
import { ref } from 'vue'
import { useElementHover } from '@vueuse/core'

const props = defineProps<{ fileName: string; index: number }>()
const emits = defineEmits<{
  view3d: [index: number]
  edit: [index: number]
  delete: [index: number]
}>()

const dom = ref()
const isHovered = useElementHover(dom)
</script>

<template>
  <div ref="dom" class="flex flex-nowrap justify-between items-center min-h-30px px-1">
    <div class="px-1">
      {{ props.fileName }}
    </div>
    <div v-if="isHovered" class="flex-(~ nowrap none) justify-between items-center gap-1">
      <div
        class="bg-indigo-2 rd-1 my-1 p-1 text-(l indigo-5) leading-1em cursor-pointer
                hover:(bg-indigo-3)
                active:(outline outline-2px outline-indigo-3)"
        @click="() => emits('view3d', props.index)"
      >
        3D
      </div>
      <div
        class="bg-sky-2 rd-1 my-1 p-1 text-(l sky-5) leading-1em cursor-pointer
                hover:(bg-sky-3)
                active:(outline outline-2px outline-sky-3)"
        @click="() => emits('edit', props.index)"
      >
        edit
      </div>
      <div
        class="bg-red-2 rd-1 my-1 p-1 text-(l red-5) leading-1em cursor-pointer
                hover:(bg-red-3)
                active:(outline outline-2px outline-red-3)"
        @click="() => emits('delete', props.index)"
      >
        del
      </div>
    </div>
  </div>
</template>

<style>
</style>
