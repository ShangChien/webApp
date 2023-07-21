<script setup lang="ts">
import { onMounted, provide, ref } from 'vue'
import { NScrollbar } from 'naive-ui'
import search from '@/components/ketcher/search.vue'
import type { pgDataItem } from '@/components/types'
import { keyMolDetail } from '@/components/types'
import gridPage from '@/components/rdkitComponent/gridPage.vue'
import molDetail from '@/components/rdkitComponent/molDetail.vue'

const mounted = ref<boolean>(false)

const result = ref<pgDataItem[]>([])

const molDetailData = { result: ref({}), searchState: ref(0) }
provide(keyMolDetail, molDetailData)
onMounted(() => {
  mounted.value = true
})
</script>

<template>
  <div class="flex flex-nowrap items-start justify-center relative box-border gap-5px h-full">
    <div class="w-60% flex-auto flex flex-col flex-nowrap items-center justify-start box-border gap-5px max-h-90vh relative">
      <div class="w-full flex-none h-full">
        <search v-model:queryResult="result" />
      </div>
      <div class="w-full flex-auto box-border">
        <grid-page :mol-list="result" :cols="6" :rows="9" class="w-full h-79vh" />
      </div>
    </div>
    <div class="w-40% flex-auto box-border b-(solid 2 indigo-100 rd-2) ">
      <NScrollbar :style="{ maxHeight: '89vh' }">
        <div v-if="molDetailData.searchState.value === 0" class="flex-(~ col) justify-center items-center box-border h-89vh">
          <div class="i-fluent-emoji-flat-ghost text-4xl text-4xl m-5" />
          <div>double-click target to get detail info</div>
        </div>
        <div v-if="molDetailData.searchState.value === 1" class="flex-(~ col) justify-center items-center box-border h-89vh">
          <div class="i-noto-hamburger animate-bounce-alt animate-count-infinite animate-duration-1.5s text-4xl" />
          <div>processing...</div>
        </div>
        <div v-if="molDetailData.searchState.value === 2" class="h-full p-1 box-border ">
          <molDetail v-bind="molDetailData.result.value" />
        </div>
        <div v-if="molDetailData.searchState.value === 3" class="flex-(~ col) justify-center items-center box-border h-89vh">
          <div class="i-fxemoji-expressionless text-4xl m-5" />
          <div>no found</div>
        </div>
        <div v-if="molDetailData.searchState.value === 4" class="flex-(~ col) justify-center items-center  box-border h-89vh">
          <div class="i-twemoji-sad-but-relieved-face text-4xl m-5" />
          <div>SORRY: Error occurred! </div>
          <div>Please repeat the operation again,</div>
          <div>Or report to developer</div>
        </div>
      </NScrollbar>
    </div>
  </div>
</template>

<style scoped>
.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.21s cubic-bezier(0.4, 0, 0.2, 1);
}

.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}
.MCD {
  animation: rotate 3s linear infinite
}
@keyframes rotate {
  0% {
      transform: rotateY(0);
  }
  25% {
      transform: rotateY(30deg);
  }
  50% {
      transform: rotateY(90deg);
  }
  75% {
      transform: rotateY(150deg);
  }
  100% {
      transform: rotateY(180deg);
  }
}
</style>
