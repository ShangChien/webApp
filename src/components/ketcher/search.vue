<script setup lang="ts">
import { computed, h, onMounted, provide, ref, toValue } from 'vue'
import type { Ref, VNodeChild } from 'vue'
import { createReusableTemplate, useDebounceFn, useElementSize } from '@vueuse/core'
import { useSortable } from '@vueuse/integrations/useSortable'
import dayjs from 'dayjs'
import { NDatePicker, NDropdown, NPopover, useMessage } from 'naive-ui'
import axios from 'axios'
import tags from '@/components/ketcher/tags.vue'
import vmodelSmiles from '@/components/ketcher/vmodelSmiles.vue'
import modalKetcher from '@/components/ketcher/modalKetcher.vue'
import { useSearchStore } from '@/stores/searchStore'
import { useMolStore } from '@/stores/molStore'
import type { condition, option, pgData, pgDataItem, recordFull, records } from '@/components/types'
import { keyStateKetcher } from '@/components/types'

enum logicType {
  And = 'and',
  Not = 'not',
  Or = 'or',
}

const ketcher = { id: ref(0), showModal: ref(false), smiles: ref('') }
provide(keyStateKetcher, ketcher)

const queryResult = defineModel<pgDataItem[] | pgData[]>()
const [DefineFilter, ReuseFilter] = createReusableTemplate<{ data: condition }>()
const [DefineTime, ReuseTime] = createReusableTemplate<{ data: {
  value: Ref<number>
  valueFormat: Ref<string>
  showDetail: Ref<boolean>
} }>()
const [DefineInputNum, ReuseInputNum] = createReusableTemplate<{ data: {
  value: Ref<number>
  placeholder: string
  unit: string
} }>()
const [DefineHistroy, ReuseHistroy] = createReusableTemplate<{ data: records[] }>()
const message = useMessage()
function info(type: any, content: string | (() => VNodeChild)) {
  message.create(
    content,
    {
      type,
      closable: true,
      duration: 3000,
      keepAliveOnHover: true,
    },
  )
}

const searchStore = useSearchStore()

// const input = ref(null)
const filter_c = ref(null)
const showFilter = ref<boolean>(true)
const canSerach = ref<boolean>(true)
const { height: h_filter_c } = useElementSize(filter_c)
const inputText = ref('')
// const showModal = ref(false)
// function clearText() {
//   inputText.value = ''
//   input.value.focus()
// }
const conditions = ref<condition[]>([])
useSortable(filter_c, conditions, {
  handle: '.handle',
  animation: 150,
})
const logic_options = [
  { label: 'and', key: logicType.And, icon: () => h('div', { class: 'i-carbon-shape-intersect bg-blue ml-1 text-2xl' }) },
  { label: 'or', key: logicType.Or, icon: () => h('div', { class: 'i-carbon-shape-unite bg-blue ml-1 text-2xl' }) },
  { label: 'not', key: logicType.Not, icon: () => h('div', { class: 'i-carbon-send-to-back bg-blue ml-1 text-2xl' }) },
]
const options = ref<option[]>([
  {
    label: '搜索方式',
    key: 'methods',
    editable: false,
    icon: () => h('div', { class: 'i-solar-card-search-bold-duotone bg-blue mr--2 text-2xl' }),
    children: [
      { label: '子结构', key: 'substructure', icon: () => h('div', { class: 'i-carbon-assembly bg-blue mr--2 text-2xl' }) },
      { label: '相似度', key: 'similarity', icon: () => h('div', { class: 'i-carbon-assembly-reference bg-blue mr--2 text-2xl' }) },
    ],
  },
  {
    label: '材料名称',
    key: 'name',
    icon: () => h('div', { class: 'i-solar-hashtag-square-bold-duotone bg-blue mr--2 text-2xl' }),
    children: [
      { label: '计算编号', key: 'name_calc', icon: () => h('div', { class: 'i-tabler-number bg-blue mr--2 text-2xl' }) },
      { label: '材料编号', key: 'name_mat', icon: () => h('div', { class: 'i-fluent-text-number-format-20-filled bg-blue mr--2 text-2xl' }) },
    ],
  },
  {
    label: '理化性质',
    key: 'properties',
    editable: true,
    icon: () => h('div', { class: 'i-solar-atom-bold-duotone bg-blue mr--2 text-2xl' }),
    children: [
      { label: '分子量', key: 'mw', icon: () => h('div', { class: 'i-ph-scales-duotone bg-blue mr--2 text-2xl' }) },
      { label: 'LUMO', key: 'lumo', icon: () => h('div', { class: 'i-ic-twotone-cloud bg-blue mr--2 text-2xl' }) },
      { label: 'HOMO', key: 'homo', icon: () => h('div', { class: 'i-ic-sharp-wb-cloudy bg-blue mr--2 text-2xl' }) },
      { label: 'Eg', key: 'eg', icon: () => h('div', { class: 'i-material-symbols-align-space-around-rounded bg-blue mr--2 text-2xl' }) },
      { label: '极化率', key: 'polar', icon: () => h('div', { class: 'i-solar-magnet-bold-duotone bg-blue mr--2 text-2xl' }) },
      { label: '电子传输', key: 'transport', icon: () => h('div', { class: 'i-solar-square-transfer-horizontal-bold-duotone bg-blue mr--2 text-2xl' }) },
    ],
  },
  {
    label: '时间范围',
    key: 'date',
    editable: false,
    icon: () => h('div', { class: 'i-solar-calendar-date-bold-duotone bg-blue mr--2 text-2xl' }),
    items: [
      { label: '开始', key: 'start' },
      { label: '结束', key: 'end' },
    ],
  },
  {
    label: '材料类型',
    key: 'types',
    editable: true,
    icon: () => h('div', { class: 'i-solar-layers-bold-duotone bg-blue mr--2 text-2xl' }),
    items: [
      { label: '蓝光客体', key: 'BD' },
      { label: '蓝光主体', key: 'BH' },
      { label: '盖层材料', key: 'CP' },
      { label: '电子阻挡', key: 'EB' },
      { label: '电子注入', key: 'EI' },
      { label: '电子传输', key: 'ET' },
      { label: '绿光客体', key: 'GD' },
      { label: '绿光主体', key: 'GH' },
      { label: '空穴注入', key: 'HI' },
      { label: '空穴传输', key: 'HT' },
      { label: '红光客体', key: 'RD' },
      { label: '红光主体', key: 'RH' },
      { label: '黄光客体', key: 'YD' },
      { label: '黄光主体', key: 'YH' },
      { label: '其他材料', key: 'etc' },
    ],
  },
  {
    label: '材料标签',
    key: 'labels',
    editable: true,
    icon: () => h('div', { class: 'i-solar-tag-price-bold-duotone bg-blue mr--2 text-2xl' }),
    items: [
      { key: 'reference', label: '对比材料' },
      { key: 'self_design', label: '自主研发' },
      { key: 'molecule', label: '分子' },
      { key: 'core', label: '主核' },
      { key: 'fragment', label: '支链' },
    ],
  },
  {
    label: '反应匹配',
    key: 'reaction',
    disabled: true,
    icon: () => h('div', { class: 'i-ion-erlenmeyer-flask bg-blue mr--2 text-2xl' }),
  },
])
const id_condition = ref<number>(0)

async function add_condition() {
  conditions.value.push({
    id: id_condition.value,
    logic: logicType.And,
    logic_icon: h('div', { class: 'i-carbon-shape-intersect bg-blue ml-1 text-2xl' }),
    label: '-select-',
    label_icon: h('div', { class: 'i-solar-alt-arrow-down-linear bg-blue text-2xl' }),
  })
  id_condition.value += 1
  showFilter.value = true
}
async function del_condition(id: number) {
  const index = await Promise.resolve(conditions.value.findIndex(el => el.id === id))
  if (index !== -1)
    conditions.value.splice(index, 1)

  else
    console.log(id, 'not found,internal error!')

  // console.log('delete',conditions.value)
}
//
const debouncedFnDay = useDebounceFn((newVal, data) => {
  const stmp = dayjs(newVal, 'YYYY.MM.DD HH:mm:ss').valueOf()
  if (!Number.isNaN(stmp))
    data.value = stmp
}, 500)
const debouncedFnText = useDebounceFn((newVal: string, data) => {
  const strList = newVal.split(/[,;\s\n]+/).map(e => e.startsWith('/') ? e : e.toUpperCase())
  data.value = strList.join(' ')
}, 500)
function onSelect(option: option, data: condition) {
  data.label = option.label as string
  data.label_icon = option.icon
  data.key = option.key
  const unit = { mw: 'M(相对分子质量)', homo: 'eV', lumo: 'eV', eg: 'eV', polar: 'C·m(Debye)', transport: 'cm²/s', similarity: '%' }
  if (option.key === 'date') {
    data.type = 'date'
    data.value = [ref(Date.now()), ref(Date.now())]
    data.showDetail = [ref(false), ref(false)]
    data.valueFormat = [
      computed<string>({
        get() {
          return dayjs(data.value[0].value).format('YYYY.MM.DD HH:mm:ss')
        },
        set(newVal) {
          debouncedFnDay(newVal, data.value[0])
        },
      }),
      computed<string>({
        get() {
          return dayjs(data.value[1].value).format('YYYY.MM.DD HH:mm:ss')
        },
        set(newVal) {
          debouncedFnDay(newVal, data.value[1])
        },
      }),
    ]
    data.component = () => h(
      'div',
      { class: 'flex items-center justify-start' },
      [h(ReuseTime, {
        data: {
          value: data.value[0],
          valueFormat: data.valueFormat[0],
          showDetail: data.showDetail[0],
        },
      }),
      h('div', { class: 'i-carbon-pan-horizontal text-gray-5' }),
      h(ReuseTime, {
        data: {
          value: data.value[1],
          valueFormat: data.valueFormat[1],
          showDetail: data.showDetail[1],
        },
      })],
    )
  }
  else if (option.key === 'labels') {
    data.type = 'tag'
    const store = useMolStore()
    const labels = computed(() => store.getAllLabels)
    data.value = []
    data.component = () => h(
      tags,
      {
        'tags': labels.value,
        'selected': data.value,
        'onUpdate:selected': (val: string[]) => { data.value = val },
      },
    )
  }
  else if (option.key === 'types') {
    data.type = 'tag'
    const labels = option.items.map(e => e.key)
    data.value = []
    data.component = () => h(
      tags,
      {
        'tags': labels,
        'selected': data.value,
        'onUpdate:selected': (val: string[]) => { data.value = val },
      },
    )
  }
  else if (option.key === 'substructure') {
    data.value = ''
    data.type = '>'
    data.component = () => h(
      'div',
      { class: 'flex items-center justify-start gap-2 rd-1 h-24px mx-1 w-full' },
      [
        h(
          'div',
          { class: 'flex-(~ basis-30% shrink-1 grow-1) items-center justify-center gap-2' },
          ['<', '=', '>'].map((item, index) => h('div',
            {
              class: [data.type === item ? 'bg-blue-3' : '',
                'bg-blue-1 text-blue-9 text-xl flex items-center justify-center h-22px w-22px box-border rd-1 cursor-pointer hover:(bg-blue-2)'],
              key: index,
              onClick: () => { data.type = item },
            },
            item,
          ),
          )),
        h(
          vmodelSmiles,
          {
            'class': 'flex-(~ basis-60% shrink-1 grow-1)',
            'conditionId': data.id,
            'smiles': data.value,
            'onUpdate:smiles': (val: string) => {
              data.value = val
              console.log(data.value, data)
            },
          },
        ),
      ],
    )
  }
  else if (option.key === 'similarity') {
    data.value = [ref<number>(), ref<number>(), ref<string>('')]
    data.type = 'interval'
    data.component = () => h(
      'div',
      { class: 'flex items-center justify-start gap-2 rd-1 h-24px mx-1 w-full' },
      [
        h(
          'div',
          { class: 'flex-(~ basis-30% shrink-1 grow-1) justify-center items-center' },
          [h(ReuseInputNum, {
            data: {
              value: data.value[0],
              placeholder: '- ∞',
              unit: '',
            },
          }),
          h('div', { class: 'i-carbon-pan-horizontal text-gray-5' }),
          h(ReuseInputNum, {
            data: {
              value: data.value[1],
              placeholder: '+ ∞',
              unit: unit[option.key],
            },
          })],
        ),
        h(
          vmodelSmiles,
          {
            'class': 'flex-(~ basis-60% shrink-1 grow-1)',
            'conditionId': data.id,
            'smiles': data.value[2].value,
            'onUpdate:smiles': (val: string) => {
              data.value[2].value = val
              console.log(data.value[2], data)
            },
          },
        ),
      ],
    )
  }
  else if (['mw', 'homo', 'lumo', 'eg', 'polar', 'transport'].includes(option.key)) {
    data.type = 'interval'
    data.value = [ref<number>(), ref<number>()]
    data.component = () => h(
      'div',
      { class: 'flex-auto flex justify-start items-center mx-1' },
      [h(ReuseInputNum, {
        data: {
          value: data.value[0],
          placeholder: '- ∞',
          unit: '',
        },
      }),
      h('div', { class: 'i-carbon-pan-horizontal text-gray-5' }),
      h(ReuseInputNum, {
        data: {
          value: data.value[1],
          placeholder: '+ ∞',
          unit: unit[option.key],
        },
      })],
    )
  }
  else if (['name_mat', 'name_calc'].includes(option.key)) {
    data.value = ref('')
    data.valueFormat = computed<string>({
      get() {
        return data.value
      },
      set(newVal) {
        debouncedFnText(newVal, data)
      },
    })
    data.component = () => h('input', {
      class: ['pl-1 rd-1 b-0 m-1 p-0 h-24px w-full outline-0 box-border bg-slate-2 hover:bg-blue-1 text-gray-8 whitespace-nowrap'],
      value: data.valueFormat,
      onInput: (e: any) => { data.valueFormat = e.target.value },
    })
  }
  else {
    // todo
  }
}

const searchField = computed<{
  key: string
  logic: string
  value: string | number | string[] | number[]
}[]>(() => {
  const out = conditions.value.reduce(
    (accumu, curr) => {
      if (curr.label !== '-select-') {
        let { key, logic, value } = curr
        if (Array.isArray(value))
          value = value.map(i => toValue(i))
        else
          value = toValue(value)

        const currValue = { key, logic, value }
        return accumu.concat(currValue)
      }
      return accumu
    },
    [],
  )
  if (inputText.value.length > 0) {
    out.push({
      key: 'smiles',
      logic: 'and',
      value: inputText.value,
    })
  }
  return out
})

function search() {
  console.log('query:', searchField.value)
  if (searchField.value.some(item => ['substructure', 'similarity'].includes(item.key))
  && !searchField.value.some(item => item.key === 'smiles')) {
    console.log('smiles is empty!')
    info('warning', '条件错误: 指定了结构搜索但是没有给出smiles式')
    return undefined
  }
  if (searchField.value.length === 0) {
    info('warning', '条件错误: 请添加搜索条件')
    return undefined
  }
  canSerach.value = false
  axios.post(
    'api/search',
    searchField.value,
  ).then(async (res: any) => {
    console.log('receve res:', res.data.data)
    queryResult.value = res.data.data
    const data4strore: recordFull | any = {
      timestamp: dayjs().unix(),
      length: queryResult.value.length,
      conditions: conditions.value,
      res: res.data.data as {
        id: number
        smiles: string
        name_calc: string
        name_mat: string
      }[],
    }
    await searchStore.addRecord(data4strore)
    canSerach.value = true
  }).catch((error) => {
    console.log(error)
    canSerach.value = true
  })
}
const font = ref<string>('font-Xl-R')
function randomPick() {
  const options = ['font-TJ-R', 'font-LXMG-R', 'font-YZ-R', 'font-Xl-R']
  const randomIndex = Math.floor(Math.random() * options.length)
  font.value = options[randomIndex]
}
const clickOut = (x: Ref<boolean>) => setTimeout(() => x.value = false, 0)
onMounted(() => {
  // check('s')
  // console.log(h_filter_c.value)
})
</script>

<template>
  <div>
    <DefineTime v-slot="{ data }">
      <NPopover
        class="p-0 m-0 rd-2" trigger="manual"
        :show="data.showDetail.value" raw
        :show-arrow="false"
        @clickoutside="clickOut(data.showDetail)"
      >
        <template #trigger>
          <div
            class="flex-auto flex items-center justify-start
              rd-1 box-border m-1 p-0 h-24px max-w-210px bg-slate-2
              hover:bg-blue-1"
            @click="data.showDetail.value = true"
          >
            <div class="flex-auto ml-1">
              <input
                v-model="data.valueFormat.value" type="text"
                class="b-0 p-0 m-0 outline-0 bg-inherit w-185px text-gray-8 text-0.9em"
              >
            </div>
            <div class="i-carbon-event-schedule bg-gray-5 mr-1 text-xl p-0 flex-none" />
          </div>
        </template>
        <NDatePicker
          v-model:value="data.value.value" panel
          type="datetime"
          value-format="yyyy.MM.dd HH:mm:ss"
          :actions="['confirm', 'now']"
          @confirm="data.showDetail.value = false"
        />
      </NPopover>
    </DefineTime>
    <DefineInputNum v-slot="{ data }">
      <input
        v-model="data.value.value"
        class="rd-1 b-0 p-0 h-24px outline-0 box-border
          w-30% max-w-210px min-w-50px text-center text-0.9em
          bg-slate-2 text-gray-8
          hover:bg-blue-1"
        :placeholder="data.placeholder"
        type="number"
        step="0.001"
      >
      <span class="text-1em text-gray-5">{{ data.unit }}</span>
    </DefineInputNum>
    <DefineFilter v-slot="{ data }">
      <div
        class="flex flex-nowrap box-border items-center justify-between bg-slate-1 gap-2 p-1 rd-2 w-full"
      >
        <div
          class="flex-none flex justify-center items-center
            h-32px w-20px box-border p-0 ml-0.5 m-0
            bg-blue-2 rd-1.5 z-1
            handle cursor-grab"
        >
          <div class="i-charm-grab-vertical flex-none bg-blue-5 hover:bg-blue-3" />
        </div>
        <div class="flex-none w-70px">
          <NDropdown
            width="trigger"
            :options="logic_options"
            :render-label="(option:any) => h('div', { class: 'w-20px' }, option.label)"
            placement="bottom-start"
            trigger="click"
            @select="(key, option) => { data.logic = key; data.logic_icon = option.icon }"
          >
            <div
              class="flex flex-nowrap justify-start items-center gap-1 transition-210 cursor-pointer
                rd-1 h-32px box-border bg-light-1 p-0 m-0 text-lg
                b-(solid 2 blue-2 rd-1) outline-(~ transparent 2px)
                hover:(b-blue-4 bg-slate-2)
                focus-within:(outline-blue-2 b-blue-4 bg-slate-2)"
              tabindex="0"
            >
              <component :is="data.logic_icon" :key="data.id" />
              <div class="p-0 m-0 text-l">
                {{ data.logic }}
              </div>
            </div>
          </NDropdown>
        </div>
        <div class="flex-none w-120px">
          <NDropdown
            width="trigger"
            :options="options"
            :render-label="(option:any) => h('div', { class: 'mr--2 ml-1' }, option.label)"
            placement="bottom-start"
            trigger="click"
            @select="(_key, option:any) => onSelect(option, data)"
          >
            <div
              class="flex flex-nowrap justify-start items-center gap-3 pl-2 transition-210 cursor-pointer
                rd-1 h-32px box-border bg-light-1 p-0 m-0 text-lg
                b-(solid 2 blue-2 rd-1) outline-(~ transparent 2px)
                hover:(b-blue-4 bg-slate-2)
                focus-within:(outline-blue-2 b-blue-4 bg-slate-2)"
              tabindex="0"
            >
              <component :is="data.label_icon" />
              <div
                :class="[data.label === '-select-' ? 'text-gray-5' : '']"
                class="p-0 m-0 text-l"
              >
                {{ data.label }}
              </div>
            </div>
          </NDropdown>
        </div>
        <div
          class="flex-auto flex items-center justify-between overflow-hidden
            h-32px box-border transition-210 bg-light-1
            b-(solid 2 blue-2 rd-1) outline-(~ transparent 2px)
            hover:b-blue-4
            active:outline-blue-2
            focus-within:(outline-blue-2 b-blue-4)"
          tabindex="0"
        >
          <component :is="data.component" />
        </div>
        <div
          class="i-solar-close-square-bold-duotone w-30px flex-none text-3xl text-blue cursor-pointer
          hover:text-blue-5"
          @click="del_condition(data.id)"
        />
      </div>
    </DefineFilter>
    <DefineHistroy v-slot="{ data }">
      <div
        class="flex flex-col flex-nowrap items-center justify-start font-LX-B
          box-border p-1 rd-2 bg-white relative min-h-400px"
      >
        <div
          class="flex-none flex flex-nowrap items-center justify-between
            rd-1 w-full min-w-400px box-border bg-slate-50"
        >
          <div class="m-1">
            {{ data.length === 0 ? 'No History' : `Total: ${data.length}` }}
          </div>
          <div
            class="bg-slate-2 rd-1 m-1 p-1 text-l leading-1em cursor-pointer
              hover:(bg-slate-3) active:(outline outline-2px outline-blue-2)"
            @click="searchStore.clearAll()"
          >
            clear all
          </div>
        </div>
        <div
          v-if="data.length === 0"
          class="flex-auto flex flex-col justify-center items-center text-2xl gap-5"
          :class="[font]"
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
            <div class="flex-basis-30 ">
              length
            </div>
            <div class="flex-basis-40 ">
              actions
            </div>
          </div>
          <div
            v-for="i in data" :key="i.id"
            class="flex-none flex flex-nowrap items-center justify-between box-border w-full p-1"
          >
            <div class="flex-basis-15">
              {{ i.id + 1 }}
            </div>
            <div class="flex-basis-60">
              {{ dayjs.unix(i.timestamp).format('YYYY-MM-DD HH:mm:ss') }}
            </div>
            <div class="flex-basis-30">
              {{ i.length }}
            </div>
            <div
              class="flex-basis-40 flex flex-nowrap items-center justify-center gap-2
                box-border rd-1 bg-lightblue-50 p-1"
            >
              <div
                class="bg-slate-2 rd-1 text-0.8em leading-4 cursor-pointer hover:b-slate-3
                  pr-1 pl-1 active:(outline outline-2px outline-blue-2)"
              >
                view
              </div>
              <div
                class="bg-red-3 rd-1 text-0.8em leading-4 cursor-pointer hover:b-red-4
                  pr-1 pl-1 active:(outline outline-2px outline-red-2)"
                @click="searchStore.rmRecordById(i.id)"
              >
                delete
              </div>
            </div>
          </div>
        </div>
      </div>
    </DefineHistroy>
    <div class="flex flex-nowrap flex-col items-center c-black gap-2 box-border w-full p-1">
      <div
        class="flex-none flex flex-nowrap justify-between items-center relative w-full
          h-36px p-0 box-border transition-210 b-(2 solid blue-4 rd-1)
          outline-(~ transparent 2px)
          focus-within:(outline-blue-2)"
        tabindex="0"
      >
        <div class="pre_line" />
        <input
          class="flex-auto text-lg outline-0 b-0 rd-1 pl-1 p-0 box-border h-full bg-slate-1 hover:(bg-blue-50)"
          placeholder="conditions"
        >
        <div
          v-if="canSerach" class="flex justify-center items-center
            box-border b-0 m--0.5 rd-r-1 p-1.5 pl-2 pr-2 bg-blue-4 cursor-pointer
            hover:(bg-blue-5) active:(outline outline-2px outline-blue-2 bg-blue-6)"
          @click="search"
        >
          <div class="i-twemoji-magnifying-glass-tilted-left text-2xl" />
        </div>
        <div
          v-else class="flex justify-center items-center
            box-border b-0 m--0.5 rd-r-1 p-1.5 pl-2 pr-2 bg-blue-4 cursor-pointer
            hover:(bg-blue-5)
            active:(outline outline-2px outline-blue-2 bg-blue-6)"
        >
          <div class="i-eos-icons-loading text-2xl" />
        </div>
      </div>
      <div
        class="flex-auto flex flex-col  w-full box-border p-0 transition-height-210"
        :style="{ height: `${h_filter_c + 40}px` }"
      >
        <div
          v-show="showFilter" ref="filter_c"
          class="flex-none flex flex-col box-border gap-2 relative"
        >
          <modalKetcher />
          <TransitionGroup name="list">
            <div v-for="item in conditions" :key="item.id" class="box-border mr-2 p-0 ">
              <ReuseFilter :data="item" />
            </div>
          </TransitionGroup>
        </div>
        <div class="flex-none flex justify-between box-border items-center bg-slate-1 rd-5 mt-2 p-1">
          <div
            class="flex-none aspect-ratio-1 bg-slate-1 rd-50% h-6 z-1 box-border
              flex justify-center items-center cursor-pointer b-(solid 2 blue-2 rd-50%)
              hover:(bg-slate-2 b-blue-3)
              active:(outline outline-2px outline-blue-2)"
          >
            <div class="i-iconamoon-sign-plus-bold text-xl text-blue" @click="add_condition()" />
          </div>
          <div
            class="flex-none flex justify-end items-center box-border mr-1 gap-1 relative"
          >
            <div
              v-show="conditions.length > 0"
              class="bg-slate-2 rd-3 box-border p-1 text-0.8em leading-3.5 cursor-pointer
                hover:(bg-slate-3) active:(outline outline-2px outline-blue-2)"
              @click="() => !(conditions.length > 0) ? null : showFilter = !showFilter"
            >
              {{ showFilter ? '折叠条件' : '展开条件' }}
            </div>
            <NPopover
              raw display-directive="if" trigger="click" :show-arrow="true" placement="right-start"
              class="p-0 rd-2 box-border"
            >
              <ReuseHistroy :data="searchStore.$state.records" />
              <template #trigger>
                <div
                  class="bg-slate-2 rd-3 box-border p-1 text-0.8em leading-3.5 cursor-pointer
                    hover:(bg-slate-3) active:(outline outline-2px outline-blue-2)"
                  @click="randomPick"
                >
                  历史记录
                </div>
              </template>
            </NPopover>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<style>
.pre_line {
  position: absolute;
  pointer-events: none;
  top:34px;
  left: 13px;
  width: 2px;
  height: v-bind("`${h_filter_c + 30}px`");
  background-color:#bfdbfe;
  z-index: 1;
}
.list-move,
.list-enter-active,
.list-leave-active {
  transition: all 0.21s ease;
}
.list-enter-from,
.list-leave-to {
  opacity: 0;
}
.list-leave-active {
  position: absolute;
}
</style>
