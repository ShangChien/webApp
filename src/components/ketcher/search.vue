<script setup lang='ts'>
import { ref, reactive, computed, onMounted } from "vue"
import { useIDBKeyval } from '@vueuse/integrations/useIDBKeyval'
import { useElementBounding,createReusableTemplate } from '@vueuse/core'
import { NButton,NModal,NCard } from "naive-ui";
import initKetcher from './initKetcher.vue';
import type { molData } from "@/components/types";
enum CharType {
  And = 'And',
  Not = 'Not',
  Or = 'Or',
}
interface option{
  key?:string;
  label?:string;
  editable?:boolean;
  value?:(number|string|boolean)[]
  items?:option[]
}
interface condition extends option{
  id:number;
  search_logic?:CharType
}


const el_condition = ref(null)
const [DefineC, ReuseC] = createReusableTemplate<{data:condition }>()

const input = ref()
const initMol: molData = reactive({qsmiles:'*~*'})
const inputText = ref("");
const showModal = ref(false);
function clearText() {
  inputText.value=''
  input.value.focus()
}
const conditions=ref<condition[]>([])
const options=ref<option[]>([
  {
    key:'methods', label:'搜索方式', editable:false,
    items:[
      {label:'子结构',key:'substructure'},
      {label:'相似度',key:'similarity'}
    ]
  },
  {label:'时间范围',key:'during',editable:false,
    items:[
      {label:'开始',key:'start'},
      {label:'结束',key:'end'}
    ]},
  {label:'理化性质',key:'properties',editable:true,
    items:[
      {label:'分子量',key:'mw'},
      {label:'HOMO',key:'homo'},
      {label:'LUMO',key:'lumo'},
      {label:'Eg',key:'eg'},
      {label:'极化率',key:'polar'},
      {label:'电子传输',key:'transport'}
  ]},
  {label:'类型',key:'types',editable:true,
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
  {key:'labels',label:'标签',editable:true,
    items:[
      {key:'reference',label:'对比材料'},
      {key:'self_design',label:'自主研发'},
      {key:'molecule',label:'分子'},
      {key:'core',label:'主核'},
      {key:'fragment',label:'支链'}
  ]},
  {key:'cal_name',label:'计算名称'},
  {key:'mat_name',label:'材料名称'}
])
const id_condition=ref<number>(0)
function add_condition() {
  conditions.value.push({id:id_condition.value})
  id_condition.value+=1
}

function search(){
  initMol.smiles = inputText.value;
  initMol.atoms = {};
  initMol.bonds = {};
  //console.log(initMol)
}
</script>
<template>
<div>
<define-c v-slot="{data}">
  <div>{{data.id}}</div>
</define-c>
<div class="flex flex-nowrap flex-col items-center text-lg c-black gap-2">
  <div class="flex-none flex flex-nowrap justify-between items-center w-60vw
    h-36px m-0.5 pr-0 rd-1 b-blue-2 box-border b-2 transition-210
    active:(outline outline-2px outline-blue-2)
    focus-within:(outline outline-2px outline-blue-2 b-blue-4)
    hover:b-blue-4" >
    <input ref="input" class='flex-auto outline-0 b-0 rd-1 pl-1 m-0 p-0 box-border h-100%'
      v-model="inputText" placeholder="smiles, smarts or Nothing"/>
    <div v-show="!!inputText" class="aspect-ratio-1 rd-50% h-5 mr-1 bg-slate-100
     flex justify-center items-center hover:(cursor-pointer bg-gray-200)"
     @click="clearText">
      <div class="i-ion-close" ></div>
    </div>
    <div @click="showModal = true"
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
  <div class="flex-none flex flex-col item-start w-60vw relative gap-2 ">
    <div v-for="item in conditions" class="b-1 b-blue-5 rd-2">
      <reuse-c :data="item"></reuse-c>
    </div>
    <div @click="add_condition()" class="i-carbon-add-alt text-2xl"></div>
  </div>
</div>
</div>
</template>
<style>
</style>