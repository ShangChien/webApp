<script setup lang="ts">
import { ref, onMounted, watch, reactive } from 'vue';
import type { molData } from '@/components/types'
import  initRDKit  from  '@/components/RDKit'
const props = defineProps<molData>()
const emit = defineEmits(['update-mol'])
const highlightMap:molData=reactive({
  smiles:props.smiles,
  atoms:[],
  bonds:[],
})
const svgitem=reactive<HTMLDivElement|any>({
  svg:{},
  rect:{},
  ellipse:[],
  path:{
    hightBonds:[],
    symble:[],
    bond:[]
  },
})
const parser = new DOMParser();
const xmlDoc=ref<Document>()
function getSvgData(svg:string){
  svgitem.svg={}
  svgitem.rect={}
  svgitem.ellipse=[]
  svgitem.path={
    hightBonds:[],
    symble:[],
    bond:[]
  }
  xmlDoc.value = parser.parseFromString(svg,"text/xml")
  //获取svgRoot内容
  let svgRoot=xmlDoc.value.getElementsByTagName('svg')[0]
  for (var attrs of ['xmlns','xmlns:rdkit','xmlns:xlink','version','baseProfile','xml:space','width','height','viewBox']){
     svgitem.svg[attrs]=svgRoot.getAttribute(attrs)
  }
   //获取svgRect内容
  let svgRect=xmlDoc.value.getElementsByTagName('rect')[0]
  for (var attrs of ['style','width','height','x','y']){
     svgitem.rect[attrs]=svgRect.getAttribute(attrs)
  }
  //获取svgEllipse内容（高亮无符号的原子）
  let svgEllipse:any=xmlDoc.value.getElementsByTagName('ellipse')
  for (var item of svgEllipse){
    svgitem.ellipse.push({ellipse:{}})
    for (var attrs of ['cx','cy','rx','ry','class']){
      svgitem.ellipse[svgitem.ellipse.length-1].ellipse[attrs]=item.getAttribute(attrs)
    }
    svgitem.ellipse[svgitem.ellipse.length-1].ellipse['style']=item.getAttribute('style').concat(';opacity: 0.6')
  }
  //获取svgPath内容
  let svgPath:any=xmlDoc.value.getElementsByTagName('path')
  for (var item of svgPath){
    //字符匹配
    if (item.getAttribute('style')==null){
      svgitem.path.symble.push({path:{}})
      for (var attrs of ['class','d','fill']){
        svgitem.path.symble[svgitem.path.symble.length-1].path[attrs]=item.getAttribute(attrs)
      }
    //化学键匹配
    } else if (item.getAttribute('style')?.split(";")[3]=='stroke-width:1.0px'){
      svgitem.path.bond.push({path:{}})
      for (var attrs of ['class','d','style']){
        svgitem.path.bond[svgitem.path.bond.length-1].path[attrs]=item.getAttribute(attrs)
      }
    //高亮化学键
    } else {
      svgitem.path.hightBonds.push({path:{}})
      for (var attrs of ['class','d','style']){
        svgitem.path.hightBonds[svgitem.path.hightBonds.length-1].path[attrs]=item.getAttribute(attrs)
      }
      svgitem.path.hightBonds[svgitem.path.hightBonds.length-1].path['style']=item.getAttribute('style').concat(';opacity: 0.2') 
    }
  }
}

function renderMol(props:molData){
    window.RDKit.prefer_coordgen(true)
    let mol= window.RDKit.get_mol(props.smiles ?? '')
    let qmol=window.RDKit.get_qmol(props.qsmiles ?? '')
    //let mDetail = JSON.parse(mol.get_substruct_match(qmol))
    let mdetailsRaw = mol.get_substruct_matches(qmol)
    let mDetail = mdetailsRaw.length > 2 ? JSON.parse(mdetailsRaw) : []
    mDetail = mDetail.reduce((acc:any, { atoms, bonds }:any) => ({atoms: [...acc.atoms, ...atoms],bonds: [...acc.bonds, ...bonds]}),{ bonds: [], atoms: [] })
    mDetail['atoms']=mDetail.atoms?.concat(props.atoms ?? []) ?? props.atoms
    mDetail['bonds']=mDetail.bonds?.concat(props.bonds ?? []) ?? props.bonds
    mDetail['addAtomIndices']=props.addAtomIndices
    mDetail['addBondIndices']=props.addBondIndices
    mDetail['lenged']=props.legend
    mDetail['width']=props.width ?? 200
    mDetail['height']=props.height ?? 200
    mDetail['highlightColour']=props.highlightColor ?? [0.624,0.675,0.902]
    mDetail['bondLineWidth']=props.bondLineWidth ?? 1
    mDetail['highlightBondWidthMultiplier']=props.highlightBondWidthMultiplier ?? 20
    mDetail['highlightRadius']=props.highlightRadius ?? 0.25
    mDetail['minFontSize']=props.minFontSize ?? 10
    mDetail['explicitMethyl']=props.explicitMethyl ?? false
    console.log(mDetail)
    mDetail=JSON.stringify(mDetail)
    let svg=mol.get_svg_with_highlights(mDetail)
    mol.delete()
    qmol.delete()
    getSvgData(svg)
    
    //console.log(bondMap.value)
    console.log(svgitem)
    //rdkitdiv.value.innerHTML=svg
}

function domClick($event:any){
  let itemList=$event.target.getAttribute('class').split(' ')[0].split('-')
  if (itemList[0]=="atom"){//如果点击atom
    //改颜色
    $event.target.style.stroke='#3fada8'
    $event.target.style.fill='#3fada8'
    $event.target.style.opacity=0.8
    //添加index到数组
    highlightMap.atoms?.push(itemList[1]/1)
    highlightMap.atoms=Array.from(new Set(highlightMap.atoms)).sort()
    emit('update-mol', highlightMap)
    //console.log(highlightMap.atoms)
  } else if(itemList[0]=="bond"){//如果点击bond
      $event.target.style.stroke='#3fada8'
      $event.target.style.fill='#3fada8'
      $event.target.style.opacity=0.4
      //添加index到数组
      highlightMap.bonds?.push(itemList[1]/1)
      highlightMap.bonds=Array.from(new Set(highlightMap.bonds)).sort()
      emit('update-mol', highlightMap)
     //console.log(highlightMap.bonds)
  }else{
    console.log(itemList[0],'error')
  }
  
}
function domDblClick($event:any){
  let itemList=$event.target.getAttribute('class').split(' ')[0].split('-')
  if (itemList[0]=="atom"){//如果点击atom
    //改颜色
    $event.target.style.stroke='#9FACE6'
    $event.target.style.fill='#9FACE6'
    $event.target.style.opacity=0.6
    //删除index
    highlightMap.atoms?.splice(highlightMap.atoms?.indexOf(itemList[1]/1),1)
    highlightMap.atoms=Array.from(new Set(highlightMap.atoms)).sort()
    emit('update-mol', highlightMap)
    //console.log(highlightMap.highlightAtoms)

  } else if(itemList[0]=="bond"){//如果点击bond
    $event.target.style.stroke='#9FACE6'
    $event.target.style.fill='#9FACE6'
    $event.target.style.opacity=0.2
    //删除index
    highlightMap.bonds?.splice(highlightMap.bonds.indexOf(itemList[1]/1),1)
    highlightMap.bonds=Array.from(new Set(highlightMap.bonds)).sort()
    emit('update-mol', highlightMap)
    //console.log(highlightMap.highlightAtoms)
  }else{
    console.log(itemList[0],'error')
  }
}

onMounted(
  async()=>{
    await initRDKit.then((res)=>{
      window.RDKit=res
    })
    renderMol(props)
  }
)

watch(
  props,
  (newVal)=>{
    renderMol(newVal)
    highlightMap.smiles=newVal.smiles
    highlightMap.atoms=[]
    highlightMap.bonds=[]
  }
)

</script>

<template class="svg">
<svg v-bind="svgitem.svg" style="width:100%; height:100%">
  <rect v-bind="svgitem.rect" />
  <path v-for="item in svgitem.path.symble" v-bind="item.path" />
  <path v-for="item in svgitem.path.bond" v-bind="item.path"  />
  <path v-for="item in svgitem.path.hightBonds" v-bind="item.path" @click="domClick($event)" @dblclick="domDblClick($event)" />
  <ellipse v-for="item in svgitem.ellipse" v-bind="item.ellipse" @click="domClick($event)" @dblclick="domDblClick($event)" />
</svg>
</template>
<style>
.svg { 
    display: inline-block;
    position: relative;
    width: 100%;
    padding-bottom: 100%; 
    vertical-align: middle; 
    overflow: hidden; 
}
</style>