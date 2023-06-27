<script setup lang="ts">
import { computed, onMounted, reactive, ref } from 'vue'
import { createReusableTemplate, useElementBounding } from '@vueuse/core'
import search from '@/components/ketcher/search.vue'
import { useMolStore } from '@/stores/molStore'
import { useSearchStore } from '@/stores/searchStore'
import type { pgDataItem } from '@/components/types'
import gridPage from '@/components/rdkitComponent/gridPage.vue'

interface dataPick {
  items: string[]
  selected: string[]
  expend: boolean
  enabled: boolean
}
const mounted = ref<boolean>(false)
const angle = ref(0)
const el_label = ref(null)
const { height: labels_h, right: labels_r } = useElementBounding(el_label)
const showExpand = computed(() => {
  return !!(el_label.value?.lastElementChild.getBoundingClientRect().right > labels_r.value
    || labels_h.value > 32)
})
const [DefineTag, ReuseTag] = createReusableTemplate<{ label: string }>()
const store = useMolStore()
const searchStore = useSearchStore()
const labels = computed(() => store.getAllLabels)
const types = ['all', 'molecule', 'ligand', 'core']
const dataLabels: dataPick = reactive({ items: labels, selected: [], expend: false, enabled: false })
const dataTypes: dataPick = reactive({ items: types, selected: ['all'], expend: false, enabled: false })
const mols4view = computed(() => {
  return store.getByTypesAndLabels(dataTypes.selected, dataLabels.selected)
})

async function select_all_or_none(data: dataPick) {
  const len = await Promise.resolve(data.items.length)
  if (len !== 0)
    data.selected = (data.selected.length !== 0) ? [] : [...data.items]

  else
    console.log('data items length', data.selected.length)
}
async function reverse(data: dataPick) {
  // 翻转icon
  angle.value = (angle.value === 0) ? 180 : 0
  const arr = await Promise.resolve(data.items)
  data.selected = arr.filter(item => !data.selected.includes(item))
}

const result = ref<pgDataItem[]>([])

onMounted(() => {
  mounted.value = true
})
</script>

<template>
  <div class="flex flex-nowrap items-center justify-center relative box-border">
    <div class="w-60% flex-auto flex flex-col flex-nowrap items-center justify-start box-border gap-5px">
      <div class="w-full flex-none">
        <search v-model="result" />
      </div>
      <div class="w-full flex-auto box-border">
        <grid-page :mol-list="result" :cols="6" :rows="9" class="w-full h-90vh" />
      </div>
    </div>
    <div class="w-40% flex-auto">
      <div v-for="(item, index) in searchStore.$state.records" :key="index" />
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
</style>
