<script setup lang="ts">
import { ref, reactive, computed,toValue,onMounted,h } from "vue";
import type { VNode,Component,Ref } from 'vue'
import { useIDBKeyval } from '@vueuse/integrations/useIDBKeyval';
import { createReusableTemplate,useElementSize,refDebounced,useDebounceFn,useFocus } from '@vueuse/core';
import { useSortable } from '@vueuse/integrations/useSortable'
import dayjs from 'dayjs'
import { NButton,NModal,NCard,NDropdown,NDatePicker,NPopover,NInputNumber } from "naive-ui"
import { DatePicker as TDatePicker } from 'tdesign-vue-next'
import tags from "@/components/ketcher/tags.vue"
import { useMolStore } from '@/stores/molStore'
import type { DropdownOption } from 'naive-ui'
import initKetcher from './initKetcher.vue'
import type { molData } from "@/components/types"
import { eachMinuteOfInterval } from "date-fns";
enum logicType {
  And = 'and',
  Not = 'not',
  Or = 'or',
}
enum filterType {
  Interval = 'Interval',
  text = 'text',
  multi_v = 'multi_v',
  radio = 'radio',
}
interface option{
  label?:string;
  editable?:boolean;
  value?:string;
  key?:string;
  icon?:VNode|Component;
  disabled?:boolean;
  items?:option[];
  children?:option[];
}
interface condition extends Omit<option,'value'>{
  id:number;
  logic:logicType;
  label?:string;
  type?:string;
  value?: string|number|string[]|number[]|Ref<number>[]|Ref<string>[]|Ref|any;
  valueFormat?:string|number|string[]|number[]|Ref<number>[]|Ref<string>[]|Ref|any;
  showDetail?:Ref<boolean>[]|Ref;
  logic_icon?:Component;
  label_icon?:VNode|Component;
  component?:Component;
}

const [DefineFilter, ReuseFilter] = createReusableTemplate<{ data:condition }>()
const [DefineTime, ReuseTime] = createReusableTemplate<{data:{
  value:Ref<number>;
  valueFormat:Ref<string>;
  showDetail:Ref<boolean>
}}>()
const [DefineInputNum, ReuseInputNum] = createReusableTemplate<{data:{
  value:Ref<number>;
  placeholder:string;
  unit:string;
}}>()
const el_condition =  ref(null)
const input        =  ref(null)
const filter_c     =  ref(null)
const { height:h_filter_c }=useElementSize(filter_c)
const initMol: molData = reactive({qsmiles:'*~*'})
const inputText = ref("")
const showModal = ref(false)
function clearText(){
  inputText.value=''
  input.value.focus()
}
const conditions=ref<condition[]>([])
useSortable(filter_c, conditions, {
  handle: '.handle',
  animation:150,
})
const logic_options=[
  {label:'and',key:logicType.And,icon:(()=>h('div',{class:'i-carbon-shape-intersect bg-blue ml-1 text-2xl'}))},
  {label:'or' ,key:logicType.Or,icon:(()=>h('div',{class:'i-carbon-shape-unite bg-blue ml-1 text-2xl'}))},
  {label:'not',key:logicType.Not,icon:(()=>h('div',{class:'i-carbon-send-to-back bg-blue ml-1 text-2xl'}))},
]
const options = ref<option[]>([
  {label:'搜索方式',key:'methods',editable:false,
    icon:(()=>h('div',{class:'i-solar-card-search-bold-duotone bg-blue mr--2 text-2xl'})),
    children:[
      {label:'子结构',key:'substructure',icon:
      (()=>h('div',{class:'i-carbon-assembly bg-blue mr--2 text-2xl'}))},
      {label:'相似度',key:'similarity',icon:(()=>h('div',{class:'i-carbon-assembly-reference bg-blue mr--2 text-2xl'}))}
    ],
  },
  {label:'材料名称',key:'name',
    icon:(()=>h('div',{class:'i-solar-hashtag-square-bold-duotone bg-blue mr--2 text-2xl'})),
    children:[
      {label:'计算编号',key:'name_cal',icon:(()=>h('div',{class:'i-tabler-number bg-blue mr--2 text-2xl'}))},
      {label:'材料编号',key:'name_mat',icon:(()=>h('div',{class:'i-fluent-text-number-format-20-filled bg-blue mr--2 text-2xl'}))},
    ]
  },
  {label:'理化性质',key:'properties',editable:true,
    icon:(()=>h('div',{class:'i-solar-atom-bold-duotone bg-blue mr--2 text-2xl'})),
    children:[
      {label:'分子量',key:'mw',icon:(()=>h('div',{class:'i-ph-scales-duotone bg-blue mr--2 text-2xl'}))},
      {label:'LUMO',key:'lumo',icon:(()=>h('div',{class:'i-ic-twotone-cloud bg-blue mr--2 text-2xl'}))},
      {label:'HOMO',key:'homo',icon:(()=>h('div',{class:'i-ic-sharp-wb-cloudy bg-blue mr--2 text-2xl'}))},
      {label:'Eg',key:'eg',icon:(()=>h('div',{class:'i-material-symbols-align-space-around-rounded bg-blue mr--2 text-2xl'}))},
      {label:'极化率',key:'polar',icon:(()=>h('div',{class:'i-solar-magnet-bold-duotone bg-blue mr--2 text-2xl'}))},
      {label:'电子传输',key:'transport',icon:(()=>h('div',{class:'i-solar-square-transfer-horizontal-bold-duotone bg-blue mr--2 text-2xl'}))}
  ]},
  {label:'时间范围',key:'date',editable:false,
  icon:(()=>h('div',{class:'i-solar-calendar-date-bold-duotone bg-blue mr--2 text-2xl'})),
    items:[
      {label:'开始',key:'start'},
      {label:'结束',key:'end'}
    ]},
  {label:'材料类型',key:'types',editable:true,
    icon:(()=>h('div',{class:'i-solar-layers-bold-duotone bg-blue mr--2 text-2xl'})),
    items:[
      {label:'蓝光客体',key:'BD',},
      {label:'蓝光主体',key:'BH',},
      {label:'盖层材料',key:'CP',},
      {label:'电子阻挡',key:'EB',},
      {label:'电子注入',key:'EI',},
      {label:'电子传输',key:'ET',},
      {label:'绿光客体',key:'GD',},
      {label:'绿光主体',key:'GH',},
      {label:'空穴注入',key:'HI',},
      {label:'空穴传输',key:'HT',},
      {label:'红光客体',key:'RD',},
      {label:'红光主体',key:'RH',},
      {label:'黄光客体',key:'YD',},
      {label:'黄光主体',key:'YH',},
      {label:'其他材料',key:'etc',},
  ]},
  {label:'材料标签',key:'labels',editable:true,
    icon:(()=>h('div',{class:'i-solar-tag-price-bold-duotone bg-blue mr--2 text-2xl'})),
    items:[
      {key:'reference',label:'对比材料'},
      {key:'self_design',label:'自主研发'},
      {key:'molecule',label:'分子'},
      {key:'core',label:'主核'},
      {key:'fragment',label:'支链'}
  ]},
  {label:'反应匹配',key:'reaction',disabled:true,
    icon:(()=>h('div',{class:'i-ion-erlenmeyer-flask bg-blue mr--2 text-2xl'})),
  }
])
const id_condition=ref<number>(0)
async function add_condition() {
  conditions.value.push({
    id:id_condition.value,
    logic:logicType.And,
    logic_icon:h('div',{class:'i-carbon-shape-intersect bg-blue ml-1 text-2xl'}),
    label:'-select-',
    label_icon:h('div',{class:'i-solar-alt-arrow-down-linear bg-blue text-2xl'}),
  })
  id_condition.value+=1
}
async function del_condition(id:number){
  let index =await Promise.resolve(conditions.value.findIndex((el)=>el.id==id))
  if (index !== -1) {
    conditions.value.splice(index, 1)
  }else{
    console.log(id,'not found,internal error!')
  }
  //console.log('delete',conditions.value)
}
//
const debouncedFnDay = useDebounceFn((newVal,data) => {
  const stmp = dayjs(newVal, 'YYYY.MM.DD HH:mm:ss').valueOf()
  if (!isNaN(stmp)) {
    data.value = stmp;
  }
}, 500)
const debouncedFnText = useDebounceFn((newVal:string,data) => {
  let strList=newVal.split(/[,;\s\n]+/) .map(e=>e.startsWith('/') ? e : e.toUpperCase())
  data.value = strList.join(' ')
}, 500)
function onSelect(option:option,data:condition) {
  data.label=option.label as string; 
  data.label_icon=option.icon;
  if (option.key=='date'){
    data.type='date'
    data.value=[ref(Date.now()),ref(Date.now())]
    data.showDetail=[ref(false),ref(false)]
    data.valueFormat=[
      computed<string>({
        get(){
          return dayjs(data.value[0].value).format('YYYY.MM.DD HH:mm:ss')
        },
        set(newVal){
          //@ts-ignore
          debouncedFn(newVal,data.value[0])
        }
      }),
      computed<string>({
        get(){
          return dayjs(data.value[1].value).format('YYYY.MM.DD HH:mm:ss')
        },
        set(newVal){
          debouncedFnDay(newVal,data.value[1])
        }
      }),
    ]
    data.component = ()=>h(
      'div',
      {class:"flex items-center justify-start"},
      [h(ReuseTime,{data:{
        value:data.value[0],
        valueFormat:data.valueFormat[0],
        showDetail:data.showDetail[0],
      }}),
      h('div',{class:'i-carbon-pan-horizontal text-gray-5'}),
      h(ReuseTime,{data:{
        value:data.value[1],
        valueFormat:data.valueFormat[1],
        showDetail:data.showDetail[1],
      }})]
    )
  }else if ('labels'==option.key) {
    data.type='tag'
    const store = useMolStore()
    const labels=computed(()=>store.getAllLabels)
    data.value=[]
    data.component = ()=>h(
      tags,
      {
        tags:labels.value,
        selected:data.value,
        'onUpdate:selected': (val:string[])=>{data.value=val},
      }
    )
  }else if ('types'==option.key) {
    data.type='tag'
    let labels=option.items.map(e=>e.key)
    data.value=[]
    data.component = ()=>h(
      tags,
      {
        tags:labels,
        selected:data.value,
        'onUpdate:selected': (val:string[])=>{data.value=val},
      }
    )
  }else if (['mw','homo','lumo','eg','polar','transport','similarity'].includes(option.key)){
    let unit={'mw':'M(相对分子质量)','homo':'eV','lumo':'eV','eg':'eV','polar':'C·m(Debye)','transport':'cm²/s','similarity':'%'}
    data.type='interval'
    data.value=[ref<number>(),ref<number>()]
    data.component= ()=>h(
      'div',
      {class:"flex-auto flex justify-start items-center"},
      [h(ReuseInputNum,{data:{
        value: data.value[0],
        placeholder:'- ∞',
        unit:'',
      }}),
      h('div',{class:'i-carbon-pan-horizontal text-gray-5'}),
      h(ReuseInputNum,{data:{
        value: data.value[1],
        placeholder:'+ ∞',
        unit:unit[option.key],
      }})]
    )
  }else if (option.key=='substructure') {
    data.value='gt'
    data.component= ()=>h(
      'div',
      {class: 'flex items-center justify-start gap-2 bg-slate-2 rd-1.2 h-24px ml-1 pr-0.5 pl-0.5'},
      [h(
        'div',
        {
          class: [toValue(data.value)=='gt' ? 'bg-blue-3 text-blue-6':'hover:bg-blue-1 text-gray-5',
            'h-22px flex items-center justify-center rd-1 pr-1 pl-1 cursor-pointer'],
          onClick: ()=>{data.value='gt'}
        },
        [h('div',{class:'i-carbon-radio-button-checked text-center'}),
        h('span',{class:'text-0.9em'},'大于检索结构')]
      ),
      h(
        'div',
        {
          class: [toValue(data.value)=='lt' ? 'bg-blue-3 text-blue-6':'hover:bg-blue-1 text-gray-5',
            'h-22px flex items-center justify-center rd-1 pr-1 pl-1 cursor-pointer'],
          onClick: ()=>{data.value='lt'}
        },
        [h('div',{class: 'i-carbon-recording-filled-alt text-center'}),
        h('span',{class: 'text-0.9em'},'小于检索结构')]
      )]
    )
  }else if (['name_mat','name_cal'].includes(option.key)){
    data.value=ref('')
    data.valueFormat = computed<string>({
      get(){
        return data.value
      },
      set(newVal){
        debouncedFnText(newVal,data)
      }
    })
    data.component= ()=>h('input',{
      class: ['pl-1 rd-1 b-0 m-1 p-0 h-24px w-full outline-0 box-border bg-slate-2 hover:bg-blue-1 text-gray-8 whitespace-nowrap'],
      value: data.valueFormat,
      'onInput': (e:any)=>{data.valueFormat=e.target.value;console.log(data.value)}
    })
  }else {
    //todo
  }
}

function search(){
  initMol.smiles = inputText.value
  initMol.atoms = {}
  initMol.bonds = {}
  console.log('search',conditions.value)
  //console.log(initMol)
}
const clickOut=(x:Ref<boolean>)=>setTimeout(()=>x.value=false,0)
onMounted(()=>{
  //check('s')
  //console.log(h_filter_c.value)
})
</script>
<template>
<div>
<define-time v-slot="{data}">
  <n-popover class="p-0 m-0 rd-2" trigger="manual" 
    :show="data.showDetail.value" raw 
    :show-arrow="false"
    @clickoutside="clickOut(data.showDetail)"	 >
    <template #trigger>
      <div class="flex-auto flex items-center justify-start
        rd-1 box-border m-1 p-0 h-24px max-w-210px bg-slate-2
        hover:bg-lightblue-1" 
        @click="data.showDetail.value=true">
        <div class="flex-auto ml-1" >
          <input type="text" class="b-0 p-0 m-0 outline-0 bg-inherit w-185px text-gray-8 text-0.9em" v-model="data.valueFormat.value">
        </div>
        <div class="i-carbon-event-schedule bg-gray-5 mr-1 text-xl p-0 flex-none" ></div>    
      </div>  
    </template>
    <n-date-picker panel type="datetime"
      value-format="yyyy.MM.dd HH:mm:ss"
      :actions="['confirm','now']"
      v-model:value="data.value.value" 
      @confirm="data.showDetail.value=false" />
  </n-popover>
</define-time>
<define-input-num v-slot="{data}">
  <input class="rd-1 b-0 m-1 p-0 h-24px outline-0 box-border
    w-30% max-w-210px min-w-50px text-center text-0.9em
    bg-slate-2 text-gray-8
    hover:bg-lightblue-1"
    :placeholder='data.placeholder'
    v-model="data.value.value"
    type="number"
    step=0.001 />
    <span class="text-1em text-gray-5">{{ data.unit }}</span>
</define-input-num>
<define-filter v-slot="{data}">
<div class="flex flex-nowrap box-border items-center justify-between bg-slate-1 
  gap-2 p-1 rd-2 w-full ">
  <div class='flex-none flex justify-center items-center 
    h-32px w-20px box-border p-0 m-0
    bg-blue-2 rd-1.5 z-1
    handle cursor-grab'>
    <div class="flex-none i-charm-grab-vertical bg-blue-5 hover:bg-blue-3"></div>
  </div>
  <div class="flex-none w-70px">
    <n-dropdown
      width="trigger"
      :options="logic_options"
      :render-label="(option:any)=>h('div',{class:'w-20px'},option.label)"
      placement="bottom-start"
      trigger="click"
      @select="(key,option)=>{data.logic=key; data.logic_icon=option.icon}">
      <div class="flex flex-nowarp justify-start items-center gap-1
        transition-210 cursor-pointer
        rd-1 h-32px b-blue-2 box-border bg-light-1 b-2 p-0 m-0 text-lg
        hover:(b-blue-4 bg-slate-2)
        focus-within:(outline outline-2px outline-blue-2 b-blue-4 bg-slate-2)"
        tabindex="0">
        <component :is="data.logic_icon" />
        <div class='p-0 m-0 text-l' >{{data.logic}}</div>
      </div>
    </n-dropdown>
  </div>
  <div class="flex-none w-120px">
    <n-dropdown
      width="trigger"
      :options="(options as any)"
      :render-label="(option:any)=>h('div',{class:'mr--2 ml-1'},option.label)"
      placement="bottom-start"
      trigger="click"
      @select="(_key,option:any)=>onSelect(option,data)">
      <div class="flex flex-nowarp justify-start pl-2 items-center gap-3
        transition-210 cursor-pointer 
        rd-1 h-32px b-blue-2 box-border bg-light-1 b-2 p-0 m-0 text-lg
        hover:(b-blue-4 bg-slate-2)
        focus-within:(outline outline-2px outline-blue-2 b-blue-4 bg-slate-2)"
        tabindex="0">
        <component :is="data.label_icon" />
        <div :class="[data.label=='-select-' ? 'text-gray-5' : '']"
          class='p-0 m-0 text-l'>{{data.label}}</div>
      </div>
    </n-dropdown>
  </div>
  <div class="flex-auto flex items-center justify-between overflow-hidden
      h-32px rd-1 b-blue-2 box-border b-2 transition-210 bg-light-1
      active:(outline outline-2px outline-blue-2)
      focus-within:(outline outline-2px outline-blue-2 b-blue-4)
      hover:b-blue-4">
      <component :is="data.component" />
  </div>
  <div class="w-30px flex-none i-solar-close-square-bold-duotone text-3xl text-blue cursor-pointer
    hover:text-blue-5"
    @click="del_condition(data.id)"></div>
</div>
</define-filter>
<div class="flex flex-nowrap flex-col items-center text-lg c-black gap-2">
  <div class="flex-none flex flex-nowrap justify-between items-center w-60vw relative 
    h-36px m-0 p-0 rd-1 b-blue-2 box-border b-2 transition-210
    active:(outline outline-2px outline-blue-2)
    focus-within:(outline outline-2px outline-blue-2 b-blue-4)
    hover:b-blue-4" >
    <div class="pre_line"></div>
    <input ref="input" class='flex-auto outline-0 b-0 rd-1 pl-1 m-0 p-0 box-border h-100% '
      v-model="inputText" placeholder="smiles, smarts or Nothing"/>
    <div v-show="!!inputText" class="aspect-ratio-1 rd-50% h-5 mr-1 bg-slate-100
     flex justify-center items-center hover:(cursor-pointer bg-gray-200)"
     @click="clearText()">
      <div class="i-ion-close" ></div>
    </div>
    <div @click="showModal=true"
      class="aspect-ratio-1 rd-50% h-7 box-border m-1 mr-2 bg-slate-100 
        flex justify-center items-center hover:(cursor-pointer bg-gray-200)">
      <div class="i-carbon-cloud-satellite text-xl bg-rose-500"></div>
    </div>
    <!--modal画板区域-->
    <n-modal v-model:show="showModal" display-directive="show">
      <n-card class="w-90vw h-96vh" 
        :bordered = "false"
        size="huge"
        role="dialog"
        aria-modal="true">
        <template #default>
          <init-ketcher
            @update-smiles="(smiles) => (inputText = smiles)"
          ></init-ketcher>
        </template>
        <template #footer>
          <div class="flex justify-center gap-2">
            <n-button @Click="showModal = false">取消</n-button>
            <n-button @Click="showModal = false">确定</n-button>
          </div>
        </template>
      </n-card>
    </n-modal>
    <div class="flex justify-center items-center b-0 m--0.5 rd-r-1 p-1.5 pl-2 pr-2 bg-blue-4
      hover:(cursor-pointer bg-blue-5)
      active:(outline outline-2px outline-blue-2 bg-blue-6)"
      @click="search" >
      <div class="i-twemoji-magnifying-glass-tilted-left text-2xl"/>
    </div>
  </div>
  <div class="flex-auto flex flex-col item-start 
    w-60vw gap-2 box-border p-0  
    transition-height-210" :style="{height:h_filter_c+40+'px'}">
    <div ref="filter_c" class="flex-none flex flex-col item-start box-border gap-2 relative">
      <TransitionGroup name="list" >
        <div v-for="item in conditions" :key="item.id" class="box-border mr-2 p-0 " >
          <reuse-filter :data="item" ></reuse-filter>
        </div>
      </TransitionGroup>
    </div>
    <div class="flex-none flex justify-start items-center">
      <div class="flex-none aspect-ratio-1 bg-slate-1 rd-50% h-6 z-1
        flex justify-center items-center b-2 b-blue-2 cursor-pointer
        hover:(bg-slate-2 b-blue-3)
        active:(outline outline-2px outline-blue-2)">
        <div ref="add_c" @click="add_condition()" 
          class="i-iconamoon-sign-plus-bold text-2xl text-blue">
        </div>
      </div>
    </div>
  </div>
</div>
</div>
</template>
<style >
.pre_line {
  position: absolute;
  pointer-events: none;
  top:34px;
  left: 11px;
  width: 2px;
  height: v-bind("(h_filter_c+20)+'px'");
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