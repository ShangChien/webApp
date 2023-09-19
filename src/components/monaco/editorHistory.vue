<script setup lang='ts'>
import { computed, ref } from 'vue'
import { NPopover } from 'naive-ui'
import { useTimeAgo } from '@vueuse/core'
import dayjs from 'dayjs'

const props = defineProps<{ records: { timestamp: number; snapshot: { text: string } }[] | any }>()
const emits = defineEmits<{
  backTo: [index:number]
  diffTo: [index:number]
  clear: []
}>()

const showPop = ref(false)
const clickOut = () => setTimeout(() => showPop.value = false, 0)
const lastStamp = computed(() => props.records[0].timestamp)
const timeAgo = useTimeAgo(lastStamp)
// const dateStr = (stamp: number) => dayjs(stamp).format('YY.MM.DD HH:mm:ss')
</script>

<template>
  <div class="flex flex-nowrap justify-center items-center gap-2">
    <div class="text-nowrap bg-slate-1 px-1 rd-1">
      {{ `modifed: ${timeAgo}` }}
    </div>
    <NPopover class="" trigger="manual" placement="bottom-end" style="padding: 0" :show="showPop" @clickoutside="clickOut()">
      <template #trigger>
        <div class="i-ic-baseline-history text-2xl bg-blue text-center cursor-pointer" @click="showPop = !showPop" />
      </template>
      <div
        class="flex flex-col flex-nowrap items-center justify-start font-LX-B
          box-border p-1 rd-2 bg-white relative min-h-400px"
      >
        <div
          class="flex-none flex flex-nowrap items-center justify-between
            rd-1 w-full min-w-400px box-border bg-slate-50"
        >
          <div class="m-1">
            {{ props.records.length === 0 ? 'No History' : `Total: ${props.records.length}` }}
          </div>
          <div
            class="bg-slate-2 rd-1 m-1 p-1 text-l leading-1em cursor-pointer
              hover:(bg-slate-3) active:(outline outline-2px outline-blue-2)"
            @click="emits('clear')"
          >
            reset
          </div>
        </div>
        <div
          v-if="props.records.length === 0"
          class="flex-auto flex flex-col justify-center items-center text-2xl gap-5"
        >
          <div class="i-fxemoji-expressionless text-5xl m-5 " />
          空空如也
        </div>
        <div
          v-else class="flex-none flex flex-col flex-nowrap items-center justify-start w-full text-center"
        >
          <div
            class="flex-none flex flex-nowrap items-center justify-between box-border w-full p-1 "
          >
            <div class="flex-basis-15 ">
              id
            </div>
            <div class="flex-basis-60 ">
              date
            </div>
            <div class="flex-basis-40 ">
              actions
            </div>
          </div>
          <div
            v-for="(i, index) in props.records" :key="i.timestamp"
            class="flex-none flex flex-nowrap items-center justify-between box-border w-full p-1"
          >
            <div class="flex-basis-15">
              {{ index + 1 }}
            </div>
            <div class="flex-basis-60">
              {{ dayjs.unix(i.timestamp).format('YYYY-MM-DD HH:mm:ss') }}
            </div>
            <div
              class="flex-basis-40 flex flex-nowrap items-center justify-center gap-2
                box-border rd-1 bg-lightblue-50 p-1"
            >
              <div
                class="bg-sky-2 rd-1 text-0.8em leading-4 cursor-pointer
                  pr-1 pl-1 active:(outline outline-2px outline-sky-2 b-sky-4)"
                @click="emits('diffTo', index)"
              >
                diff view
              </div>
              <div
                class="bg-red-3 rd-1 text-0.8em leading-4 cursor-pointer
                  pr-1 pl-1 active:(outline outline-2px outline-red-2 b-red-4)"
                @click="emits('backTo', index)"
              >
                back to
              </div>
            </div>
          </div>
        </div>
      </div>
    </NPopover>
  </div>
</template>

<style>
</style>
