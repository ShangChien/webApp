<script setup lang="ts">
import { ref, onMounted, watch, reactive } from 'vue';
import type { molData } from '@/components/types'
import  initRDKitModule from "@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js"
const props = defineProps<molData>()
//const rdkitdiv=ref<HTMLDivElement|any>()

//const selectedData=reactive({
//  smiles:'',
//  atoms:[],
//  bonds:[],
//})
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
  //获取svgPath内容
  let svgPath:any=xmlDoc.value.getElementsByTagName('path')
  for (var item of svgPath){
    if (item.getAttribute('class')==null){
      svgitem.path.hightBonds.push({path:{}})
      for (var attrs of ['d','style','fill']){
        svgitem.path.hightBonds[svgitem.path.hightBonds.length-1].path[attrs]=item.getAttribute(attrs)
      }
    } else if (item.getAttribute('style')==null){
      svgitem.path.symble.push({path:{}})
      for (var attrs of ['class','d','fill']){
        svgitem.path.symble[svgitem.path.symble.length-1].path[attrs]=item.getAttribute(attrs)
      }
    }else {
      svgitem.path.bond.push({path:{}})
      for (var attrs of ['class','d','style']){
        svgitem.path.bond[svgitem.path.bond.length-1].path[attrs]=item.getAttribute(attrs)
      }
    }
  }
  //获取svgEllipse内容（高亮无符号的原子）
  let svgEllipse:any=xmlDoc.value.getElementsByTagName('ellipse')
  for (var item of svgEllipse){
    svgitem.ellipse.push({ellipse:{}})
    for (var attrs of ['cx','cy','rx','ry','style']){
      svgitem.ellipse[svgitem.ellipse.length-1].ellipse[attrs]=item.getAttribute(attrs)
    }
  }
}

function renderMol(props:molData){
  initRDKitModule().then((instance:any)=>{
  const RDKit= instance
  let mol= RDKit.get_mol(props.smiles ?? '')
  let qmol=RDKit.get_qmol(props.qsmiles ?? '')
  console.log(props.smiles)
  let mDetail= JSON.parse(mol.get_substruct_match(qmol))
  mDetail['atoms']=mDetail.atoms?.concat(props.atoms ?? []) ?? props.atoms
  mDetail['bonds']=mDetail.bonds?.concat(props.bonds ?? []) ?? props.bonds
  mDetail['addAtomIndices']=props.addAtomIndices
  mDetail['addBondIndices']=props.addBondIndices
  mDetail['lenged']=props.legend
  mDetail['width']=props.width ?? 100
  mDetail['height']=props.height ?? 100
  mDetail['highlightColour']=props.highlightColor ?? [0.624,0.675,0.902]
  mDetail['bondLineWidth']=props.bondLineWidth ?? 1
  mDetail['highlightBondWidthMultiplier']=props.highlightBondWidthMultiplier ?? 15
  mDetail['highlightRadius']=props.highlightRadius ?? 0.3
  mDetail['minFontSize']=props.minFontSize ?? 10
  mDetail['explicitMethyl']=props.explicitMethyl ?? false
  mDetail=JSON.stringify(mDetail)
  let svg=mol.get_svg_with_highlights(mDetail)
  getSvgData(svg)



  //rdkitdiv.value.innerHTML=svg
  mol.delete()
  qmol.delete()
  })
}

function domClick($event:any){
  $event.target.style.stroke='#3fada8'
  $event.target.style.fill='#3fada8'
}

onMounted(()=>{
 renderMol(props)
})

watch(
  props,
  (newVal,oldVal)=>{
    renderMol(newVal)
  },
  {deep:true}
)

</script>

<template>
<!--
  <div ref="rdkitdiv"></div>
-->
<svg v-bind="svgitem.svg" >
  <rect v-bind="svgitem.rect" />
  <ellipse v-for="item in svgitem.ellipse" v-bind="item.ellipse" @click="domClick($event)" />
  <path v-for="item in svgitem.path.hightBonds" v-bind="item.path" @click="domClick($event)" />
  <path v-for="item in svgitem.path.symble" v-bind="item.path" />
  <path v-for="item in svgitem.path.bond" v-bind="item.path"  />
</svg>
</template>