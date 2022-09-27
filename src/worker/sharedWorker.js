// import rdkit into the worker, and export the renderMol function
import initRDKitModule from "@rdkit/rdkit/dist/RDKit_minimal.js";
import { optimize } from 'svgo/lib/svgo.js';

//atoms和bonds的索引处理
const preHandleProps = (props)=>{
  const concatIndex = (list1)=>{ 
    let outList= []
    for(let i in list1){
      outList = outList.concat(list1[i])
    }
    return Array.from(new Set(outList)) 
  }
  props.atoms = concatIndex(props.atoms)
  props.bonds = concatIndex(props.bonds)
  return props
}
//css高亮颜色处理
const Color = (n) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,60%,1)'
//rdkit的renderMol函数
const renderMol = (props)=>{
  let atomsIndex = Object.values(props.atoms??{}).reduce((pre,cur)=>pre.concat(cur),[])
  let bondsIndex = Object.values(props.bonds??{}).reduce((pre,cur)=>pre.concat(cur),[])
  let mol = self.RDKit.get_mol(props.smiles ?? "");
  let qmol = self.RDKit.get_qmol(props.qsmiles ?? "");
  let mdetailsRaw = mol.get_substruct_matches(qmol);
  let mDetail = mdetailsRaw.length > 2 ? JSON.parse(mdetailsRaw) : [];
  mDetail = mDetail.reduce(
    (acc, { atoms, bonds }) => ({
      atoms: [...acc.atoms, ...atoms],
      bonds: [...acc.bonds, ...bonds],
    }),
    { bonds: [], atoms: [] }
  );
  mDetail["atoms"] = mDetail.atoms?.concat(atomsIndex ?? []) ?? atomsIndex;
  mDetail["bonds"] = mDetail.bonds?.concat(atomsIndex ?? []) ?? bondsIndex;
  mDetail["addAtomIndices"] = props.addAtomIndices;
  mDetail["addBondIndices"] = props.addBondIndices;
  mDetail["lenged"] = props.legend;
  mDetail["width"] = props.width ?? 200;
  mDetail["height"] = props.height ?? 200;
  mDetail["highlightColour"] =
    props.highlightColor ?? [208, 78, 214].map((i) => i / 255);
  mDetail["bondLineWidth"] = props.bondLineWidth ?? 1;
  mDetail["highlightBondWidthMultiplier"] =
    props.highlightBondWidthMultiplier ?? 15;
  mDetail["highlightRadius"] = props.highlightRadius ?? 0.25;
  mDetail["minFontSize"] = props.minFontSize ?? 10;
  mDetail["explicitMethyl"] = props.explicitMethyl ?? false;
  mDetail = JSON.stringify(mDetail);
  out = mol.get_svg_with_highlights(mDetail);
  //console.log(out)
  mDetail = null;
  mdetailsRaw = null;
  qmol.delete();
  mol.delete();
  return out
}
//初始化高亮原子和化学键
const initHighlight=(svg)=>{
  //console.log(svg)
  let strList =svg.split(/>\n</g)
  //console.log(strList)
  let out=strList.map((str)=>{
    if (str.match(/ellipse/)){
      let atomIndex = str.match(/atom-\d+/)[0].split('-')[1]
      let colorType = Object.keys(props?.atoms??{})
                      .find((type)=>props.atoms[type].includes(+atomIndex));
      if (colorType!==undefined){
        return str.replace(/#[A-z,0-9]{6}/g,Color(+colorType))
      }else{
        //console.log('Not matched atom svg element:',str)
        return str
      }
    } else if (str.match(/path[\s,\S]*4.8px/)) {
      let bondIndex = str.match(/bond-\d+/)[0].split('-')[1]
      let colorType = Object.keys(props?.bonds??{})
                      .find((type)=>props.bonds[type].includes(+bondIndex));
      if (colorType!==undefined){
        return str.replace(/#[A-z,0-9]{6}/g,Color(+colorType))
      }else{
        //console.log('Not matched bond svg element:',str)
        return str
      }
    } else {
      return str
    }  
  })
  return out.join('>\n<')
}
//内联到css背景中
const cssBgSvg=(svgI)=>{
  return 'url("data:image/svg+xml;utf8,'+optimize(svgI, {}).data
    .replace('<svg', (~svgI.indexOf('xmlns') ? '<svg' : '<svg xmlns="http://www.w3.org/2000/svg"'))
    .replace(/"/g, '\'')
    .replace(/%/g, '%25')
    .replace(/#/g, '%23')
    .replace(/{/g, '%7B')
    .replace(/}/g, '%7D')
    .replace(/</g, '%3C')
    .replace(/>/g, '%3E')
    .replace(/atom/g,'a')
    .replace(/bond/g,'b')+'") no-repeat center center'
}

let i = 0 //counter for the number of times the worker is connected
let out = null //the output(svg image) of the renderMol function
var props
self.onconnect = async function(e) {
  //initialize the rdkit
  if (!self.RDKit) {
    await initRDKitModule().then(res=>{
      self.RDKit = res	
      self.RDKit.prefer_coordgen(true);
    }).catch(err=>{
      console.log(err)
    })
  }
  var port = e.ports[0];
  console.log("the ",i++,"th time connected")
  port.onmessage = function(e) {
    props = JSON.parse(e.data)
    if (props) {
      //console.log('props',props)
      out = initHighlight(renderMol(props))
      //console.log('out',out)
      if (props.css){
        port.postMessage(cssBgSvg(out));
      }else{
        port.postMessage(cssBgSvg(out));
      }
      //console.log(cssBgSvg(out))
		  out=null;
    }
  }
  port.onmessageerror = function(e) {
    console.log("error", e);
  }
}
self.onerror = function(e) {
  console.log("error", e);
}






