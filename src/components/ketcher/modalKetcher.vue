<script setup lang='ts'>
import { inject, ref, toValue, watch } from 'vue'
import type { Ref } from 'vue'
import { NButton, NCard, NModal } from 'naive-ui'
import initKetcher from './initKetcher.vue'
import { keyStateKetcher } from '@/components/types'

const { showModal, smiles } = inject<{
  showModal: Ref<boolean>
  smiles: Ref<string>
}>(keyStateKetcher, {
  showModal: ref(false),
  smiles: ref(''),
})
const ketcherRef = ref(null)
watch(showModal,
  (newVal, _oldVal) => {
    if (newVal && toValue(ketcherRef.value?.mounted)) {
      toValue(ketcherRef).sendMessage()
    }
  }, {
    flush: 'post',
  },
)
function confirm() {
  smiles.value = ''
  toValue(ketcherRef).getMessage()
  const unwatch = watch(
    () => smiles.value,
    () => {
      showModal.value = false
      unwatch()
    },
  )
}
</script>

<template>
  <NModal v-model:show="showModal" display-directive="show">
    <NCard
      class="w-80vw h-80vh relative"
      :bordered="false"
      size="huge"
      role="dialog"
      aria-modal="true"
    >
      <template #default>
        <init-ketcher
          ref="ketcherRef"
          v-model:smiles="smiles"
          class="w-full h-full"
        />
      </template>
      <template #footer>
        <div class="flex justify-center gap-2">
          <NButton @click="showModal = false">
            取消
          </NButton>
          <NButton @click="confirm">
            确定
          </NButton>
        </div>
      </template>
    </NCard>
  </NModal>
</template>

<style>
</style>
