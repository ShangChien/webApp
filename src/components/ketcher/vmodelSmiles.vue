<script setup lang='ts'>
import { inject, ref, watch } from 'vue'
import { computedEager } from '@vueuse/core'
import type { Ref } from 'vue'
import { keyStateKetcher } from '@/components/types'

const props = defineProps<{ conditionId: number }>()
const smiles = defineModel<string>('smiles')
const { id, showModal, smiles: modalSmiles } = inject<{
  id: Ref<number>
  showModal: Ref<boolean>
  smiles: Ref<string>
}>(keyStateKetcher, {
  id: ref(0),
  showModal: ref(false),
  smiles: ref(''),
})
function setShowModal() {
  id.value = props.conditionId
  showModal.value = true
  modalSmiles.value = smiles.value
}

const sameId = computedEager(() => id.value === props.conditionId)
const diffSmiles = computedEager(() => modalSmiles.value !== smiles.value)
const needUpdate = computedEager(() => {
  if (sameId.value && showModal.value) {
    if (diffSmiles.value)
      return true
    else
      return false
  }
  else {
    return false
  }
})
watch(needUpdate, () => {
  smiles.value = modalSmiles.value
})

const inputRef = ref(null)
</script>

<template>
  <div
    class="flex-auto flex flex-nowrap justify-between items-center gap-2
    h-24px m-0 p-0 box-border w-full"
  >
    <input
      ref="inputRef"
      v-model="smiles"
      class="flex-auto text-lg outline-0 b-0 rd-1 pl-1 m-0 p-0 box-border h-full bg-slate-2
      hover:bg-blue-1"
      placeholder="smiles"
    >
    <div
      v-show="!!smiles"
      class="aspect-ratio-1 rd-50% h-5 bg-slate-200 flex justify-center items-center
      hover:(cursor-pointer bg-gray-300)"
      @click="() => { smiles = '';inputRef.focus() }"
    >
      <div class="i-ion-close" />
    </div>
    <div
      class="aspect-ratio-1 rd-50% h-6 box-border bg-slate-100
      flex justify-center items-center hover:(cursor-pointer bg-gray-200)"
      @click="setShowModal()"
    >
      <div class="i-carbon-cloud-satellite text-xl bg-rose-500" />
    </div>
  </div>
</template>

<style>
</style>
