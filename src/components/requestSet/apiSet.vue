<script setup lang='ts'>
import { computed, inject, ref, watch } from 'vue'
import { refDebounced } from '@vueuse/core'
import { NInput, NPopover } from 'naive-ui'
import type { Ref } from 'vue'

const apiPrefix: Ref<string> = inject('apiPrefix', ref('https://192.168.2.233:5055'))
const apiCustom = ref<string>('192.168.2.233:5056')
const debouncedUrl = refDebounced(apiCustom, 2000)
const apiPreset = '192.168.2.233:5056'

const isTyping = computed(() => apiCustom.value !== debouncedUrl.value)
const isChecking = ref(false)
const isOKapiUrl = ref(true)

watch(
  debouncedUrl,
  (val) => {
    isChecking.value = true
    isOKapiUrl.value = false

    const requestOptions: RequestInit = {
      method: 'GET',
      headers: {
        'Accept': '*/*',
        'Accept-Language': 'zh-CN,zh;q=0.9,en;q=0.8',
        'Sec-Fetch-Dest': 'empty',
        'Sec-Fetch-Mode': 'cors',
        'Sec-Fetch-Site': 'cross-site',
      },
      referrerPolicy: 'strict-origin-when-cross-origin',
      body: null,
      mode: 'cors',
      credentials: 'omit',
    }

    fetch(`https://${val}/`, requestOptions)
      .then(response => response.json())
      .then((data) => {
        if (data?.detail === 'Not Found') {
          isOKapiUrl.value = true
          apiPrefix.value = `https://${val}`
          console.log('预检请求成功')
        } else {
          isOKapiUrl.value = false
          console.log('预检请求失败')
        }
        isChecking.value = false
      })
      .catch((error) => {
        isChecking.value = false
        console.log('发生错误', error)
      })
  },
)
</script>

<template>
  <div>
    <NPopover
      raw display-directive="if" trigger="click" :show-arrow="true" placement="right-end"
      class="p-2 rd-2 bg-white box-border"
    >
      <div class="p-2 w-400px flex flex-col gap-2">
        <div class="text-lg text-nowrap flex flex-nowrap justify-between items-center">
          <div>Backend request setting: </div>
          <div
            class="text-sm rd-2 p-1 bg-blue-1 hover:(bg-blue-2 cursor-pointer)"
            @click="() => { apiPrefix = 'https://192.168.2.233:5055'; apiCustom = '192.168.2.233:5055' }"
          >
            Preset
          </div>
        </div>
        <NInput v-model:value="apiCustom" size="small" round :placeholder="apiPreset" clearable>
          <template #prefix>
            https://
          </template>
          <template #suffix>
            <div>
              <div v-if="isChecking || isTyping" class="i-svg-spinners-pulse-3 text-yellow text-lg" />
              <div v-else-if="isOKapiUrl" class="i-ep-success-filled text-green text-lg" />
              <div v-else class="i-ep-warning-filled text-red text-lg" />
            </div>
          </template>
        </NInput>
      </div>
      <template #trigger>
        <div
          class="rd-1 ma p-1 cursor-pointer h-full w-full box-border flex
          hover:(bg-slate-100 cursor-pointer duration-210 ease-in-out)"
        >
          <div class="i-solar-settings-bold-duotone text-3xl ma text-slate-5" />
        </div>
      </template>
    </NPopover>
  </div>
</template>

<style>
</style>
