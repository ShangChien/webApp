<script setup lang='ts'>
import { defineAsyncComponent, h } from 'vue'
import type { VNode } from 'vue'
import type { pgDataItem } from '@/components/types'

const props = defineProps<pgDataItem>()
const nglViewer = defineAsyncComponent({
  loader: () => import('@/components/nglViewer.vue'),
  loadingComponent: () => h(
    'div',
    { class: 'flex-(~ col) justify-center items-center box-border h-full' },
    [h('div', { class: 'i-noto-hamburger animate-bounce-alt animate-count-infinite animate-duration-1.5s text-4xl' }), 'processing...'],
  ),
  delay: 2000, // +600*Math.random(),
  errorComponent: () => h(
    'div',
    { class: 'flex-(~ col) justify-center items-center box-border h-full' },
    h('div', { class: 'i-twemoji-sad-but-relieved-face text-4xl m-5' }),
  ),
  timeout: 3000,
  suspensible: false,
})
const sortOrder = ['structures', 'id', 'timestamp', 'name_calc', 'name_mat', 'smiles', 'functional', 'basisset', 'charge', 'spinmultiplicity', 'numberbasisfunc', 'elections', 'ispcm', 'isspin', 'stationarytype', 'extrainfo', 'homo', 'lumo', 'eg', 'energy', 'sfp', 'bfp', 'polar', 'dirpath', 'corrections', 'spectrum', 'uv', 'uv_low', 'iso', 'aniso', 'polarizability', 'ground', 'udft', 'tdft', 'polar', 'zeroPoint', 'energy', 'enthalpy', 'gibbsFreeEnergy', 's0_s1', 's0_t1', 'oscillator']

function isLeaf(obj: any): boolean {
  return typeof obj !== 'object' || obj === null || Array.isArray(obj)
}
function leafDom(v: number[] | string | number | null) {
  if (['number', 'string'].includes(typeof v) || (v === null)) {
    if ((typeof v === 'string') && v.startsWith('\n     RDKit          3D\n\n')) {
      return h(nglViewer, { data: v })
    }
    return v
  }
  if (Array.isArray(v)) {
    return h('div', { class: '' }, v.join(', '))
  }
}
// C1C=CC=C(C2C(C3C=CC=CC=3)=C(C3C=CC=CC=3)C(C3C=CC=CC=3)=CC=2)C=1
function convert2Dom(k: string, v: any): VNode | (() => VNode) {
  return h(
    'div',
    { class: 'flex items-center justify-start gap-0.5 bg-indigo-1 box-border w-full rd-1' },
    [h('div', { class: 'flex-none rd-1 p-1' }, `${k}: `),
      h('div', { class: 'flex-auto break-all rd-1 bg-fuchsia-1 m-0.3 p-1 box-border' }, leafDom(v) || h('div', { class: 'i-fluent-emoji-flat-japanese-free-of-charge-button text-xl' }))],
  )
}
function renderItems(obj: any): VNode | (() => VNode) {
  const sortedKeys = Object.keys(obj).sort((a, b) => {
    const indexA = sortOrder.indexOf(a)
    const indexB = sortOrder.indexOf(b)
    // 如果 a 在 sortOrder 中的索引小于 b，则 a 应该排在 b 前面
    if (indexA < indexB) {
      return -1
    }
    // 如果 a 在 sortOrder 中的索引大于 b，则 a 应该排在 b 后面
    if (indexA > indexB) {
      return 1
    }
    // 如果 a 和 b 在 sortOrder 中的索引相同，则它们的顺序不变
    return 0
  })
  return h(
    'div',
    { class: 'flex-(~ col nowrap) items-start justify-start gap-1 w-full' },
    sortedKeys.map((k) => {
      if (isLeaf(obj[k])) {
        return h(convert2Dom(k, obj[k]))
      } else {
        return h(
          'div',
          { class: 'flex items-center justify-start gap-1 b-(solid 0.4 blue-50 rd-1) bg-blue-50 box-border w-full' },
          [h('div', {}, k), h(renderItems(obj[k]))],
        )
      }
    }),
  )
}
</script>

<template>
  <component :is="renderItems(props)" />
</template>

<style>
</style>
