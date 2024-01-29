<script setup lang='ts'>
import { computed, inject, ref, watch } from 'vue'
import { useFetch } from '@vueuse/core'
import type { Ref } from 'vue'
import type { spectrumData } from './dataCheckStore'

const props = defineProps<{
  data4check: { name: string; nm: number[]; intensity: number[] }[]
}>()
const data4refernce = defineModel<spectrumData>('data4refernce')

const apiPrefix: Ref<string> = inject('apiPrefix')

const getFigNameList = useFetch(
  `${apiPrefix.value}/get_fig_name_list`,
).get().json<string[]>()
const figNameList = computed<string[]>(() => getFigNameList.data.value ?? [])

const nameinDB = ref({ name: '' })
const getFigByName = useFetch(
  `${apiPrefix.value}/get_uv_data`,
  { immediate: false },
).post(nameinDB).json<spectrumData>()
watch(() => nameinDB.value.name, async () => {
  await getFigByName.execute()
  data4refernce.value = getFigByName.data.value
})

const nameinDB4del = ref({ name: '' })
const delFigByName = useFetch(
  `${apiPrefix.value}/del_uv_data`,
  { immediate: false },
).post(nameinDB4del)

const checkFig = useFetch(
  `${apiPrefix.value}/check_uv_data`,
  { immediate: false },
).post(props.data4check)
</script>

<template>
  <div>
    <p v-for="(item, index) in figNameList" :key="index" @click="() => { nameinDB4del.name = item; delFigByName.execute(); }">
      {{ item }}
    </p>
    <input v-model="nameinDB.name" type="text" name="nameinDB">
    <button @click="checkFig.execute()">check</button>
    <p />
  </div>
</template>

<style>
</style>
