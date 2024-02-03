<script setup lang='ts'>
import { computed, inject, ref, watch } from 'vue'
import { useFetch } from '@vueuse/core'
import type { Ref } from 'vue'
import type { Spectrum, SpectrumFromDB } from './dataCheckStore'

const props = defineProps<{ data4check: Spectrum[] }>()
const data4ref = defineModel<Spectrum[]>('data4ref')

const apiPrefix: Ref<string> = inject('apiPrefix')

const getNames = useFetch(
  `${apiPrefix.value}/spectrum/get_names`,
).get().json<{ data: string[] }>()
const figNames = computed<string[]>(() => getNames.data.value.data ?? [])

const name4get = ref('')
const getByName = useFetch(
  `${apiPrefix.value}/spectrum/get`,
  { immediate: false },
).post({ name: name4get }).json<SpectrumFromDB>()
watch(() => name4get.value, async () => {
  await getByName.execute()
  // data4ref.value = getByName.data.value
})

const nameinDB4del = ref({ name: '' })
const delByName = useFetch(
  `${apiPrefix.value}/spectrum/del`,
  { immediate: false },
).post(nameinDB4del)

const check = useFetch(
  `${apiPrefix.value}/spectrum/check`,
  { immediate: false },
).post(props.data4check)
</script>

<template>
  <div>
    <p v-for="(item, index) in figNames" :key="index" @click="() => { nameinDB4del.name = item; delByName.execute(); }">
      {{ item }}
    </p>
    <input v-model="name4get" type="text" name="nameinDB">
    <button @click="check.execute()">check</button>
    <p />
  </div>
</template>

<style>
</style>
