// import rdkit into the worker, and export the renderMol function
import initRDKitModule from "@rdkit/rdkit/dist/RDKit_minimal.js";
import { DOMParser, XMLSerializer } from "@xmldom/xmldom"
import { optimize } from 'svgo/lib/svgo.js';
//initialize the rdkit
initRDKitModule().then(res=>{
	self.RDKit = res	
  self.RDKit.prefer_coordgen(true);
}).catch(err=>{
	console.log(err)
})

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
const Color = (n) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,70%,1)'
//rdkit的renderMol函数
const renderMol=(props)=>{
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
  mDetail["atoms"] = mDetail.atoms?.concat(props.atoms ?? []) ?? atomsIndex;
  mDetail["bonds"] = mDetail.bonds?.concat(props.bonds ?? []) ?? bondsIndex;
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
  mDetail["highlightRadius"] = props.highlightRadius ?? 0.3;
  mDetail["minFontSize"] = props.minFontSize ?? 10;
  mDetail["explicitMethyl"] = props.explicitMethyl ?? false;
  mDetail = JSON.stringify(mDetail);
  out = mol.get_svg_with_highlights(mDetail);
  mDetail = null;
  mdetailsRaw = null;
  qmol.delete();
  qmol = null;
  mol.delete();
  mol = null;
  return out
}
//初始化高亮原子和化学键
const initHighlight=(svg)=>{
  var parser = new DOMParser();
  let xmlSvg = parser.parseFromString(svg, "text/xml");
  console.log(props)
  Object.keys(props?.atoms).forEach((key) => {
    props.atoms?.[key].forEach((n)=>{
      let ellipse = xmlSvg.getElementsByTagName('ellipse').findIndex((i)=>i.getAttribute('class').split('-')[1]==n.toString())
      ellipse.style.fill=Color(key)
      ellipse.style.stroke=Color(key)
      console.log(ellipse)
    })
  })
  Object.keys(props?.bonds).forEach((key) => {
    props.atoms?.[key].forEach((n)=>{
      let pathHighlight = xmlSvg.getElementsByTagName('path').findIndex((i)=>i.getAttribute('class').split('-')[1]==n.toString() && i.style.stroke=='#9FACE6')
      pathHighlight.style.stroke=Color(key)
      console.log(pathHighlight)
    })
  })
  var serializer = new XMLSerializer()
  return serializer.serializeToString(xmlSvg)
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
    .replace(/bond/g,'b')+'")'
}

let i = 0 //counter for the number of times the worker is connected
let out = null //the output(svg image) of the renderMol function
var props
self.onconnect = function(e) {
  var port = e.ports[0];
  console.log("the ",i++,"th time connected")
  port.onmessage = function(e) {
    props = JSON.parse(e.data)
    if (props) {
      out = initHighlight(renderMol(preHandleProps(props)))
      console.log('out',out)
      if (props.css){
        port.postMessage(cssBgSvg(out));
      }else{
        port.postMessage(out);
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






