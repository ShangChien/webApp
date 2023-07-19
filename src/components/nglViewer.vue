<script setup lang="ts">
import { onMounted, ref } from 'vue'
import * as NGL from 'ngl/dist/ngl.js'

const props = defineProps<{ data: string }>()
const viewport = ref(null)
function initViewer(viewport: any) {
  const stage = new NGL.Stage(viewport, { backgroundColor: 'white' })
  const stringBlob = new Blob([props.data], { type: 'text/plain' })
  stage.loadFile(stringBlob, { ext: 'sdf' }).then((o: any) => {
    o.addRepresentation('licorice')
    o.stage.setSpin(true)
    o.autoView()
    // 悬停事件
  })
  stage.signals.hovered.add((pickingProxy) => {
    if (pickingProxy) {
      stage.setSpin(false)
    } else {
      stage.setSpin(true)
    }
  })
}
onMounted(() => {
  initViewer(viewport.value)
})
</script>

<template>
  <div ref="viewport" class="aspect-ratio-square w-full b-(solid 2 indigo-100 rd-1)" />
</template>

<style>
</style>
