<script setup lang='ts'>
import { computed, inject, ref, watch } from 'vue'
import { NButton, NCheckbox, NCheckboxGroup, NCollapse, NCollapseItem, NInputNumber, NModal, NScrollbar, NSwitch, NTabPane, NTabs, NTag, useMessage } from 'naive-ui'
import { useElementSize, useFetch } from '@vueuse/core'
import type { Ref } from 'vue'
import type { Spectrum, SpectrumFromDB } from './dataCheckStore'
import { data2Spectrum } from './dataCheckStore'

const props = defineProps<{ data4check: Spectrum[] }>()
const data4ref = defineModel<Spectrum[]>('data4ref')

const message = useMessage()

const domRef = ref()
const domSize = useElementSize(domRef)

const apiPrefix: Ref<string> = inject('apiPrefix')

const getNames = useFetch(
  `${apiPrefix.value}/data_check/spectrum/get_names`,
  { immediate: true },
).get().json<{ data: string[] }>()
const figNames = computed<string[]>(() => getNames.data.value?.data ?? [])

const name4get = ref({ name: '' })
const getByName = useFetch(
  `${apiPrefix.value}/data_check/spectrum/get`,
  { immediate: false },
).post(name4get).json<{ data: SpectrumFromDB }>()
watch(() => name4get.value.name, async () => {
  if (name4get.value.name !== '') {
    await getByName.execute()
    console.log(getByName.data.value)
    data4ref.value.push(...data2Spectrum(getByName.data.value.data))
  }
})

const canDelete = ref(false)
const showAuthModal = ref(false)
const modalInput = ref<number>()
const authCheck = computed<boolean>(() => {
  const now = new Date()
  const hours = now.getHours().toString().padStart(2, '0')
  const minutes = now.getMinutes().toString().padStart(2, '0')
  const passwd = `${hours}${minutes}`.split('').map(i => (+i + 1) % 10).join('')
  return (+passwd - modalInput.value) < 5
})
function handleModeChange() {
  if (canDelete.value) {
    showAuthModal.value = false
    canDelete.value = false
  } else {
    showAuthModal.value = true
    // positive-click的回调中验证auth, 再赋值canDelete.value
  }
}
function handleAuthPositiveClick() {
  console.log(`auth check:${authCheck.value}${modalInput.value}`)
  if (authCheck.value) {
    message.success('auth success!', { duration: 4000 })
    canDelete.value = true
    showAuthModal.value = false
    modalInput.value = null
  } else {
    message.error('auth failed!', { duration: 5000 })
    modalInput.value = null
  }
}

const nameinDB4del = ref({ name: '' })
const delByName = useFetch(
  `${apiPrefix.value}/data_check/spectrum/del`,
  { immediate: false },
).post(nameinDB4del).json<{ data: string }>()
async function del(name: string) {
  nameinDB4del.value.name = name
  await delByName.execute().then(() => {
    console.log(`delete action in server db: ${delByName.data.value.data}`)
    getNames.data.value = { data: getNames.data.value?.data.filter((item: string) => item !== name) }
    console.log(`delete action in browser: ${name}`)
    message.success(`${name} delete success!`, { duration: 5000 })
  })
}

const showOverWriteModal = ref(false)
const similaritydata = ref<{ name: string; similarity: number }[]>()
const data4checkRef = computed<SpectrumFromDB>(() => {
  return {
    name: props.data4check[0].name.split('.')[0],
    raw_arr: props.data4check.map(i => i.intensity),
  }
})
const checkSimilarity = useFetch(
  `${apiPrefix.value}/data_check/spectrum/check`,
  { immediate: false },
).post(data4checkRef).json<{ data: { name: string; similarity: number }[] }>()
function handleOverWritePositiveClick() {
  checkSimilarity.execute().then(() => {
    similaritydata.value = checkSimilarity.data.value.data.sort((a, b) => b.similarity - a.similarity)
    message.success(`data-${data4checkRef.value.name} success overwrited; success check similarity!`, { duration: 4000 })
  })
}
async function check() {
  await getNames.execute()
  if (figNames.value.includes(data4checkRef.value.name)) {
    showOverWriteModal.value = true
    return
  }
  checkSimilarity.execute().then(() => {
    similaritydata.value = checkSimilarity.data.value.data.sort((a, b) => b.similarity - a.similarity)
    message.success(`data-${data4checkRef.value.name} success added; success check similarity!`, { duration: 4000 })
  })
}

const checkedNames = ref<string[]>([])
function handleCheckedGroup(value: string[], meta: { actionType: 'check' | 'uncheck'; value: string }) {
  checkedNames.value = value // 组件上无法使用v-model, 需要这里直接更新选中的数组
  if (meta.actionType === 'check') {
    console.log(value, meta)
    name4get.value.name = meta.value
  } else {
    name4get.value.name = ''
    data4ref.value = data4ref.value.filter((item: Spectrum) => item.name !== meta.value)
  }
}

function result4vis(name: string) {
  if (checkedNames.value.includes(name)) {
    // 删除该name
    name4get.value.name = ''
    checkedNames.value = checkedNames.value.filter((item: string) => item !== name)
    data4ref.value = data4ref.value.filter((item: Spectrum) => item.name !== name)
  } else {
    name4get.value.name = name
    checkedNames.value.push(name)
  }
}

function railStyle({ focused, checked }) {
  const style: any = {}
  if (checked) {
    style.background = '#d03050'
    if (focused) {
      style.boxShadow = '0 0 0 2px #d0305040'
    }
  } else {
    style.background = '#2080f0'
    if (focused) {
      style.boxShadow = '0 0 0 2px #2080f040'
    }
  }
  return style
}
</script>

<template>
  <div ref="domRef" class="m-1 h-full box-border">
    <NTabs type="segment" animated>
      <NTabPane name="chap1" tab="数据管理" display-directive="show:lazy">
        <div class="flex flex-col justify-start items-center max-h-full">
          <div class="flex-none flex justify-between items-center px-1 w-full box-border">
            <NButton type="info" secondary size="small" @click="getNames.execute()">
              <div class="i-solar-refresh-circle-bold-duotone text-xl mr-1 ml--1 inline-block" />刷新获取所有表名
            </NButton>
            <div>
              mode:
              <NSwitch :value="canDelete" :rail-style="railStyle" @update:value="handleModeChange">
                <template #checked>
                  delete
                </template>
                <template #unchecked>
                  view
                </template>
              </NSwitch>
              <NModal
                v-model:show="showAuthModal"
                preset="dialog"
                positive-text="确认"
                negative-text="取消"
                @positive-click="handleAuthPositiveClick"
                @negative-click="showAuthModal = false"
              >
                <template #icon>
                  <div class="i-solar-danger-bold-duotone text-2xl text-red text-center" />
                </template>
                <template #header>
                  切换为删除模式，危险操作，请输入授权码：
                </template>
                <template #default>
                  <NInputNumber
                    v-model:value="modalInput"
                    :show-button="false"
                    placeholder="auth code"
                  />
                </template>
              </NModal>
            </div>
          </div>
          <div class="flex-auto h-70% relative w-full">
            <NCollapse :default-expanded-names="['1']">
              <NCollapseItem name="1" display-directive="show">
                <template #header>
                  <div class="flex justify-start items-center w-full">
                    <div class="flex-none m-1">当前选取:</div>
                    <div class="flex flex-auto overflow-hidden justify-start items-center gap-1 m-1 w-50px">
                      <NTag v-for="item, key in checkedNames" :key="key" type="success" size="small" round>
                        {{ item }}
                      </NTag>
                    </div>
                  </div>
                </template>
                <template #default>
                  <NScrollbar :style="`max-height:${domSize.height.value - 120}px`">
                    <NCheckboxGroup :value="checkedNames" @update:value="handleCheckedGroup">
                      <div class="flex justify-start items-center gap-1 flex-wrap p-1">
                        <div v-for="(item, index) in figNames" :key="index" class="bg-gray-1 px-1.5 py-0.8 m-1 rd-1 flex justify-center items-center hover:bg-slate-2 cursor-pointer relative">
                          <NCheckbox :value="item" :label="item" />
                          <div v-if="canDelete" class="i-ic-baseline-cancel text-xl text-red-4 hover:text-red-5 absolute top--6px right--6px" @click="del(item)" />
                        </div>
                      </div>
                    </NCheckboxGroup>
                  </NScrollbar>
                </template>
              </NCollapseItem>
            </NCollapse>
          </div>
        </div>
      </NTabPane>
      <NTabPane name="chap2" tab="数据对比" display-directive="show:lazy">
        <div class="flex flex-col justify-start items-center gap-1 w-full box-border">
          <div class="flex-none flex justify-between items-center px-1 w-full box-border">
            <NButton type="info" secondary size="small" @click="check">检查相似度</NButton>
            <NModal
              v-model:show="showOverWriteModal"
              preset="dialog"
              positive-text="确认覆盖"
              negative-text="取消"
              @positive-click="handleOverWritePositiveClick"
              @negative-click="showOverWriteModal = false"
            >
              <template #icon>
                <div class="i-solar-danger-bold-duotone text-2xl text-red-4 text-center" />
              </template>
              <template #header>
                数据库中已存在待检查的同名文件，<br>请确认是否使用新数据覆盖旧数据!
              </template>
              <template #default>
                tip:如需更改当前文件名，请在浏览器中删除，并修改文件名字重新上传至浏览器
              </template>
            </NModal>
          </div>
          <div class="flex flex-col justify-start items-center gap-1 box-border w-full p-1">
            <div class="flex justify-between items-center gap-1 m-1 rd-1 w-full box-border bg-indigo-1">
              <div class="flex-none w-40% m-1 box-border">名称</div>
              <div class="flex-auto m-1 box-border">相似度</div>
            </div>
            <NScrollbar :style="`max-height:${domSize.height.value - 100}px`">
              <div class="flex flex-col justify-start items-center gap-1 p-1 w-full box-border">
                <div
                  v-for="i in similaritydata" :key="i.name"
                  class="flex justify-between items-center gap-1 px-1 rd-1 w-full box-border cursor-pointer  hover:bg-slate-2"
                  :class="[checkedNames.includes(i.name) ? 'bg-green-2' : 'bg-slate-50']"
                  @click="result4vis(i.name)"
                >
                  <div class="flex-none w-40%">{{ i.name }}</div>
                  <div class="flex-auto">{{ i.similarity }}</div>
                </div>
              </div>
            </NScrollbar>
          </div>
        </div>
      </NTabPane>
    </NTabs>
  </div>
</template>

<style scoped>
:deep(.n-collapse-item__content-inner) {
  padding-top: 0 !important/* 修改折叠内容的 padding */
}
</style>
