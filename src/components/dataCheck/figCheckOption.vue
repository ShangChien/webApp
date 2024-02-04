<script setup lang='ts'>
import { computed, inject, ref, watch } from 'vue'
import { NButton, NCheckbox, NCheckboxGroup, NScrollbar, NTag } from 'naive-ui'
import { useFetch } from '@vueuse/core'
import type { Ref } from 'vue'
import type { Spectrum, SpectrumFromDB } from './dataCheckStore'

const props = defineProps<{ data4check: Spectrum[] }>()
const data4ref = defineModel<Spectrum[]>('data4ref')

const apiPrefix: Ref<string> = inject('apiPrefix')

const getNames = useFetch(
  `${apiPrefix.value}/data_check/spectrum/get_names`,
  { immediate: false },
).get().json<{ data: string[] }>()
const figNames = computed<string[]>(() => getNames.data.value?.data ?? [])

const name4get = ref('')
const getByName = useFetch(
  `${apiPrefix.value}/data_check/spectrum/get`,
  { immediate: false },
).post({ name: name4get }).json<SpectrumFromDB>()
watch(() => name4get.value, async () => {
  await getByName.execute()
  // data4ref.value = getByName.data.value
})

const nameinDB4del = ref({ name: '' })
const delByName = useFetch(
  `${apiPrefix.value}/data_check/spectrum/del`,
  { immediate: false },
).post(nameinDB4del)

const checkSimilarity = useFetch(
  `${apiPrefix.value}/data_check/spectrum/check`,
  { immediate: false },
).post(props.data4check)

const checkedNames = ref<string[]>([])
function onChecked(_val, meta) {
  if (meta.actionType === 'check') {
    console.log(meta)
    nameinDB4del.value.name = meta.value
    delByName.execute()
  } else {
    console.log(meta)
  }
  // nameinDB4del.value.name = meta.value
  // delByName.execute()
}
</script>

<template>
  <div>
    <div class="flex justify-around m-1">
      <NButton type="info" secondary size="small" @click="getNames.execute()">获取所有表名</NButton>
      <NButton type="info" secondary size="small" @click="checkSimilarity.execute()">检查相似度</NButton>
    </div>
    <div class="flex justify-start items-center gap-1 m-1">
      <span>当前选取:</span>
      <NTag v-for="item, key in checkedNames" :key="key" type="success" size="small" round>{{ item }}</NTag>
    </div>
    <NScrollbar style="max-height: 700px" class="m-1">
      <NCheckboxGroup v-model:value="checkedNames" @update:value="onChecked">
        <NCheckbox v-for="(item, index) in figNames" :key="index" :value="item" :label="item" />
      </NCheckboxGroup>
    </NScrollbar>
  </div>
</template>

<style>
</style>
