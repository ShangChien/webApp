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

const data = { ...props }

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
// C1C=C(C2C=C(C3C=CC=CC=3)C(C3C=CC(C4C=CC=CC=4)=CC=3)=C(C3C=CC=CC=3)C=2)C=CC=1
function convert2Dom(k: string, v: any): VNode | (() => VNode) {
  return h(
    'div',
    { class: 'flex items-center justify-start gap-0.5 bg-indigo-1 box-border w-full' },
    [h('div', { class: 'flex-none rd-1 p-1' }, `${k}: `),
      h('div', { class: 'flex-auto break-all rd-1 bg-fuchsia-1 m-0.3 p-1' }, leafDom(v) || h('div', { class: 'i-fluent-emoji-flat-japanese-free-of-charge-button text-xl' }))],
  )
}
function renderItems(obj: any): VNode | (() => VNode) {
  return h(
    'div',
    { class: 'flex-(~ col nowrap) items-start justify-start gap-1 w-full' },
    Object.keys(obj).map((k) => {
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
  <component :is="renderItems(data)" />
</template>

<style>
</style>
