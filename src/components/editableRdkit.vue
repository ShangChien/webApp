<script setup lang="ts">
import { ref, onMounted, watch, reactive } from 'vue';
import type { molData } from '@/components/types'
import  initRDKitModule from "@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js"
const props = defineProps<molData>()
const emit = defineEmits(['update-mol'])
const highlightMap:molData=reactive({
  smiles:props.smiles,
  atoms:[],
  bonds:[],
})
const bondLengthRange=ref([0,0])
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
const bondMap=ref<any>({})
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
    for (var attrs of ['cx','cy','rx','ry']){
      svgitem.ellipse[svgitem.ellipse.length-1].ellipse[attrs]=item.getAttribute(attrs)
    }
    svgitem.ellipse[svgitem.ellipse.length-1].ellipse['style']=item.getAttribute('style').concat(';opacity: 0.6')
  }
  //获取svgPath内容
  let svgPath:any=xmlDoc.value.getElementsByTagName('path')
  for (var item of svgPath){
    if (item.getAttribute('class')==null){
      svgitem.path.hightBonds.push({path:{}})
      for (var attrs of ['d','style','fill']){
        svgitem.path.hightBonds[svgitem.path.hightBonds.length-1].path[attrs]=item.getAttribute(attrs)
      }
      svgitem.path.hightBonds[svgitem.path.hightBonds.length-1].path['style']=item.getAttribute('style').concat(';opacity: 0.2')
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
}
//bondsite和class的映射关系
function getBondList(){
  let bondList:object|any={}
  //添加所有的bond
  for (var item of svgitem.path.bond){
    let key=item.path['class']//item.path['class'].split(' ')[1] 
    let start=item.path['d'].split(' ')[1].split(',')
    start=start.map((x:any)=>x/1)
    let end=item.path['d'].split(' ')[3].split(',')
    end=end.map((x:any)=>x/1)
    //console.log(start,end)
    if (bondList[key]==null){//单键
      bondList[key]=[start,end]
    } else if (bondList[key][1][0]==start[0] && bondList[key][1][1]==start[1]){
      bondList[key][1]=end
    } else {//否则为双键
      if (bondList[key][2]==null){
        bondList[key].push(start,end)
      } else {
        if (bondList[key][3][0]==start[0] && bondList[key][3][1]==start[1]){
          bondList[key][3]=end
        } 
      }
    }
  }
  //整理合并多值bondList
  for (var key in bondList){
    if (bondList[key].length==4){
      let lengthBond1=Math.sqrt(Math.pow(bondList[key][0][0]-bondList[key][1][0],2)+Math.pow(bondList[key][0][1]-bondList[key][1][1],2))
      let lengthBond2=Math.sqrt(Math.pow(bondList[key][2][0]-bondList[key][3][0],2)+Math.pow(bondList[key][2][1]-bondList[key][3][1],2))
      if  ((lengthBond1-lengthBond2)>(lengthBond1+lengthBond2)*0.05){
        bondLengthRange.value[1]=lengthBond1
        bondLengthRange.value[0]=lengthBond2
        bondList[key]=[[bondList[key][0][0]/1,bondList[key][0][1]/1],[bondList[key][1][0]/1,bondList[key][1][1]/1]]
      } else if((lengthBond2-lengthBond1)>(lengthBond1+lengthBond2)*0.05){
        bondLengthRange.value[0]=lengthBond1
        bondLengthRange.value[1]=lengthBond2
        bondList[key]=[[bondList[key][2][0]/1,bondList[key][2][1]/1],[bondList[key][3][0]/1,bondList[key][3][1]/1]]
      }else{
        bondList[key]=[[(bondList[key][0][0]/1+bondList[key][2][0]/1)/2,(bondList[key][0][1]/1+bondList[key][2][1]/1)/2],
                       [(bondList[key][1][0]/1+bondList[key][3][0]/1)/2,(bondList[key][1][1]/1+bondList[key][3][1]/1)/2]]
      }
      //console.log(lengthBond1,lengthBond2,bondList[key],bondLengthRange.value)
    }
  }
  console.log(bondList)
  return bondList
}

function renderMol(props:molData){
 initRDKitModule().then((instance:any)=>{
    const RDKit= instance
    let mol= RDKit.get_mol(props.smiles ?? '')
    let qmol=RDKit.get_qmol(props.qsmiles ?? '')
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
    mDetail['highlightBondWidthMultiplier']=props.highlightBondWidthMultiplier ?? 30
    mDetail['highlightRadius']=props.highlightRadius ?? 0.4
    mDetail['minFontSize']=props.minFontSize ?? 10
    mDetail['explicitMethyl']=props.explicitMethyl ?? false
    mDetail=JSON.stringify(mDetail)
    let svg=mol.get_svg_with_highlights(mDetail)
    mol.delete()
    qmol.delete()
    getSvgData(svg)
    bondMap.value=getBondList()
    //console.log(bondMap.value)
    //console.log(svgitem)
    //rdkitdiv.value.innerHTML=svg
  })
}
//位点匹配
function siteMatch(x:number,y:number): any {
  for (var item in bondMap.value){
    let distance1 = Math.sqrt(Math.pow(bondMap.value[item][0][0]-x,2)+Math.pow(bondMap.value[item][0][1]-y,2))
    let distance2 = Math.sqrt(Math.pow(bondMap.value[item][1][0]-x,2)+Math.pow(bondMap.value[item][1][1]-y,2))
    if (distance2>distance1){
     if (distance1<(bondLengthRange.value[0]*0.2)){
       //console.log(item.split(' ')[1].split('-')[1],[x,y])
      return item.split(' ')[1].split('-')[1]
      } 
    } else { 
      if (distance2<(bondLengthRange.value[1]*0.2)){
         //console.log(item.split(' ')[2].split('-')[1],[x,y])
        return item.split(' ')[2].split('-')[1]
      } 
    }
  }
  for (var item in bondMap.value){
    let distance1 = Math.sqrt(Math.pow(bondMap.value[item][0][0]-x,2)+Math.pow(bondMap.value[item][0][1]-y,2))
    let distance2 = Math.sqrt(Math.pow(bondMap.value[item][1][0]-x,2)+Math.pow(bondMap.value[item][1][1]-y,2))
      if (distance2<(bondLengthRange.value[1]*0.4)){
         //console.log(item.split(' ')[2].split('-')[1],[x,y])
        return item.split(' ')[2].split('-')[1]
      } 
  }
}
//bond位点匹配
function bondSiteMatch(x1:number,y1:number,x2:number,y2:number): any {
  for (var item in bondMap.value){
    let distance1 = Math.sqrt(Math.pow(bondMap.value[item][0][0]-x1,2)+Math.pow(bondMap.value[item][0][1]-y1,2))
    let distance2 = Math.sqrt(Math.pow(bondMap.value[item][1][0]-x2,2)+Math.pow(bondMap.value[item][1][1]-y2,2))
    if ((distance2+distance1)<bondLengthRange.value[0]*0.3){
      return item.split(' ')[0].split('-')[1]
      } 
  }
  for (var item in bondMap.value){
    let distance1 = Math.sqrt(Math.pow(bondMap.value[item][0][0]-x2,2)+Math.pow(bondMap.value[item][0][1]-y2,2))
    let distance2 = Math.sqrt(Math.pow(bondMap.value[item][1][0]-x1,2)+Math.pow(bondMap.value[item][1][1]-y1,2))
    if ((distance2+distance1)<bondLengthRange.value[0]*0.3){
      return item.split(' ')[0].split('-')[1]
      } 
  }
}

function domClick($event:any){
  let atomIndex=siteMatch($event.target.getAttribute('cx')/1,$event.target.getAttribute('cy')/1)/1
  if (!isNaN(atomIndex)){//如果点击atom
    //改颜色
    $event.target.style.stroke='#3fada8'
    $event.target.style.fill='#3fada8'
    $event.target.style.opacity=0.8
    //添加index到数组
    highlightMap.atoms?.push(atomIndex)
    highlightMap.atoms=Array.from(new Set(highlightMap.atoms)).sort()
    emit('update-mol', highlightMap)
    //console.log(highlightMap.highlightAtoms)

  } else {//如果点击bond
    let line =$event.target.getAttribute('d')?.split(' ')
    //console.log(line)
    let bondIndex=bondSiteMatch(line[1].split(',')[0]/1,line[1].split(',')[1]/1,
                                line[3].split(',')[0]/1,line[3].split(',')[1]/1)/1
    if (!isNaN(bondIndex)){
      //改颜色
      $event.target.style.stroke='#3fada8'
      $event.target.style.fill='#3fada8'
      $event.target.style.opacity=0.4
      //添加index到数组
      highlightMap.bonds?.push(bondIndex)
      highlightMap.bonds=Array.from(new Set(highlightMap.bonds)).sort()
      emit('update-mol', highlightMap)
      //console.log(highlightMap.highlightBonds)
    }
  }
  
}
function domDblClick($event:any){
  let atomIndex=siteMatch($event.target.getAttribute('cx')/1,$event.target.getAttribute('cy')/1)/1
  if (!isNaN(atomIndex)){
    if (highlightMap.atoms?.indexOf(atomIndex)!=-1){
      //恢复颜色
      $event.target.style.stroke='#9FACE6'
      $event.target.style.fill='#9FACE6'
      $event.target.style.opacity=0.6
      //删除index
      highlightMap.atoms?.splice(highlightMap.atoms?.indexOf(atomIndex),1)
      highlightMap.atoms=Array.from(new Set(highlightMap.atoms)).sort()
      emit('update-mol', highlightMap)
      //console.log(highlightMap.highlightAtoms)
    }
  }else {
    let line =$event.target.getAttribute('d')?.split(' ')
    let bondIndex=bondSiteMatch(line[1].split(',')[0]/1,line[1].split(',')[1]/1,
                                  line[3].split(',')[0]/1,line[3].split(',')[1]/1)/1
    if (highlightMap.bonds?.indexOf(bondIndex)!=-1){
      //恢复颜色
      $event.target.style.stroke='#9FACE6'
      $event.target.style.fill='#9FACE6'
      $event.target.style.opacity=0.2
      //删除index
      highlightMap.bonds?.splice(highlightMap.bonds.indexOf(bondIndex),1)
      highlightMap.atoms=Array.from(new Set(highlightMap.atoms)).sort()
      emit('update-mol', highlightMap)
      //console.log(highlightMap.highlightBonds)
    }
  }
}

onMounted(()=>{

 renderMol(props)
 
})

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

<template>
<svg v-bind="svgitem.svg" >
  <rect v-bind="svgitem.rect" />
  <path v-for="item in svgitem.path.symble" v-bind="item.path" />
  <path v-for="item in svgitem.path.bond" v-bind="item.path"  />
  <path v-for="item in svgitem.path.hightBonds" v-bind="item.path" @click="domClick($event)" @dblclick="domDblClick($event)" />
  <ellipse v-for="item in svgitem.ellipse" v-bind="item.ellipse" @click="domClick($event)" @dblclick="domDblClick($event)" />
</svg>
</template>