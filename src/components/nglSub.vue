<script setup lang="ts">
import { onMounted, ref } from 'vue'
import * as NGL from 'ngl/dist/ngl.js'

const viewport = ref(null)

function initViewer(viewport: any) {
  const stage = new NGL.Stage(viewport, { backgroundColor: 'white' })
  stage.loadFile('/moleculeData/mol1.sdf').then((o: any) => {
    o.addRepresentation('licorice')
    o.stage.setSpin(true)
    o.autoView()
    console.log(o)
    // 悬停事件
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
