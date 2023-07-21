<script setup lang="ts">
import { onMounted, onUnmounted, ref, watch } from 'vue'
import * as NGL from 'ngl/dist/ngl.js'
import { useElementHover } from '@vueuse/core'

const props = defineProps<{ data: string }>()
const stringBlob = new Blob([props.data], { type: 'text/plain' })
const viewport = ref(null)
const unwatch = ref(null)
const isHovered = useElementHover(viewport)
function initViewer(viewport: any) {
  const stage = new NGL.Stage(viewport, { backgroundColor: 'white' })
  return stage
}
onMounted(() => {
  const viewer = initViewer(viewport.value)
  viewer.loadFile(stringBlob, { ext: 'sdf' }).then((o: any) => {
    o.addRepresentation('licorice')
    o.stage.setSpin(true)
    o.autoView()
    // 悬停事件
  })
  unwatch.value = watch(isHovered,
    () => {
      if (isHovered.value) {
        viewer.setSpin(false)
      } else {
        viewer.setSpin(true)
      }
    },
  )
})
onUnmounted(() => {
  unwatch.value()
})
</script>

<template>
  <div ref="viewport" class="aspect-ratio-square w-full m-0 p-0 box-border rd-1" />
</template>

<style>
</style>
