<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import * as NGL from 'ngl/dist/ngl.js'
import { useElementHover, useElementSize } from '@vueuse/core'

const props = defineProps<{ data: string }>()
const stringBlob = new Blob([props.data], { type: 'text/plain' })
const viewerWindow = ref(null)
const isHovered = useElementHover(viewerWindow)
const viewport = ref(null)
const viewportSize = useElementSize(viewerWindow)
// const stage = ref(null)
const spin = ref(false)
const autoView = ref(false)
const fullScreen = ref(false)

onMounted(() => {
  const stage = new NGL.Stage(viewport.value, { backgroundColor: 'white' })
  stage.loadFile(stringBlob, { ext: 'sdf' })
    .then((o: any) => {
      o.addRepresentation('licorice')
      // o.stage.setSpin(true)
      o.autoView()
      // 悬停事件
    })
  watch([viewportSize.width, viewportSize.height],
    () => {
      stage.handleResize()
    },
  )
  watch(spin, val => stage.setSpin(val))
  watch(autoView, () => stage.autoView(500))
  watch(fullScreen, () => stage.toggleFullscreen())
})
</script>

<template>
  <div ref="viewerWindow" class="relative h-full">
    <div class="aspect-ratio-square w-full m-0 p-0 box-border rd-1" />
    <div ref="viewport" class="absolute top-0 left-0 h-full w-full m-0 p-0 box-border rd-1" />
    <div v-show="isHovered" class="absolute bottom-0 left-0 flex flex-nowrap gap-1 items-end justify-end box-border w-full animate-fade-in animate-duration-0.3s">
      <div
        class="bg-slate-3 rd-1 m-1 p-1 text-lg leading-1em cursor-pointer i-material-symbols-3d-rotation-outline-rounded
          hover:(bg-slate-5)
          active:(outline outline-2px outline-blue-4)"
        @click="spin = !spin"
      />
      <div
        class="bg-slate-3 rd-1 m-1 p-1 text-lg leading-1em cursor-pointer i-material-symbols-view-in-ar-outline-sharp
          hover:(bg-slate-5)
          active:(outline outline-2px outline-blue-4)"
        @click="autoView = !autoView"
      />
      <div
        class="bg-slate-3 rd-1 m-1 p-1 text-lg leading-1em cursor-pointer i-ic-sharp-open-in-full
          hover:(bg-slate-5)
          active:(outline outline-2px outline-blue-4)"
        @click="fullScreen = !fullScreen"
      />
    </div>
  </div>
</template>

<style>
</style>
