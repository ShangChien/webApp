/* eslint-disable */
// import rdkit into the worker, and export the renderMol function
import initRDKitModule from '@rdkit/rdkit/dist/RDKit_minimal.js'
import { optimize } from 'svgo/lib/svgo.js'

let i = 0 // counter for the number of times the worker is connected
let out = null // the output(svg image) of the renderMol function
let props
// css高亮颜色处理
const Color = n => `hsla(${(n + 8.6) * 36},90%,60%,1)`
// rdkit的renderMol函数
function renderMol(props) {
  const atomsIndex = Object.values(props.atoms ?? {}).reduce((pre, cur) => pre.concat(cur), [])
  const bondsIndex = Object.values(props.bonds ?? {}).reduce((pre, cur) => pre.concat(cur), [])
  const mol = self.RDKit.get_mol(props.smiles ?? '')
  const qmol = self.RDKit.get_qmol(props.qsmiles ?? '')
  let mdetailsRaw = mol.get_substruct_matches(qmol)
  let mDetail = mdetailsRaw.length > 2 ? JSON.parse(mdetailsRaw) : []
  mDetail = mDetail.reduce(
    (acc, { atoms, bonds }) => ({
      atoms: [...acc.atoms, ...atoms],
      bonds: [...acc.bonds, ...bonds],
    }),
    { bonds: [], atoms: [] },
  )
  mDetail.atoms = mDetail.atoms?.concat(atomsIndex ?? []) ?? atomsIndex
  mDetail.bonds = mDetail.bonds?.concat(bondsIndex ?? []) ?? bondsIndex
  mDetail.addAtomIndices = props.addAtomIndices
  mDetail.addBondIndices = props.addBondIndices
  mDetail.lenged = props.legend
  mDetail.width = props.width ?? 200
  mDetail.height = props.height ?? 200
  mDetail.highlightColour
    = props.highlightColor ?? [208, 78, 214].map(i => i / 255)
  mDetail.bondLineWidth = props.bondLineWidth ?? 1
  mDetail.highlightBondWidthMultiplier
    = props.highlightBondWidthMultiplier ?? 15
  mDetail.highlightRadius = props.highlightRadius ?? 0.25
  mDetail.minFontSize = props.minFontSize ?? 10
  mDetail.explicitMethyl = props.explicitMethyl ?? false
  mDetail = JSON.stringify(mDetail)
  out = mol.get_svg_with_highlights(mDetail)
  console.log(out)
  mDetail = null
  mdetailsRaw = null
  qmol.delete()
  mol.delete()
  return out
}
// 初始化高亮原子和化学键
function initHighlight(svg) {
  // console.log(svg)
  const strList = svg.split(/>\n</g)
  let bondsLen = Object.keys(props?.bonds ?? {}).reduce((sum, cur) => sum + props.bonds[cur].length, 0)
  let atomsLen = Object.keys(props?.atoms ?? {}).reduce((sum, cur) => sum + props.atoms[cur].length, 0)
  // console.log(strList)
  const out = strList.map((str) => {
    if ((str.search(/ellipse/) !== -1) && (atomsLen >= 1)) {
      atomsLen -= 1
      const atomIndex = str.match(/atom-\d+/)[0].split('-')[1]
      const colorType = Object.keys(props?.atoms ?? {})
        .find(type => props.atoms[type].includes(+atomIndex))
      if (colorType !== undefined) {
        return str.replace(/#[A-z,0-9]{6}/g, Color(+colorType))
      }
      else {
        // console.log('Not matched atom svg element:',str)
        return str
      }
    }
    else if ((str.search(/path[\s,\S]*bond-\d+\satom-\d/) !== -1) && (bondsLen >= 1)) {
      bondsLen -= 1
      const bondIndex = str.match(/bond-\d+/)[0].split('-')[1]
      const colorType = Object.keys(props?.bonds ?? {})
        .find(type => props.bonds[type].includes(+bondIndex))
      if (colorType !== undefined) {
        return str.replace(/#[A-z,0-9]{6}/g, Color(+colorType))
      }
      else {
        // console.log('Not matched bond svg element:',str)
        return str
      }
    }
    else {
      return str
    }
  })
  return out.join('>\n<')
}
// 内联到css背景中
function cssBgSvg(svgI) {
  return `url("data:image/svg+xml;utf8,${optimize(svgI, {}).data
    .replace('<svg', (~svgI.indexOf('xmlns') ? '<svg' : '<svg xmlns="http://www.w3.org/2000/svg"'))
    .replace(/"/g, '\'')
    .replace(/%/g, '%25')
    .replace(/#/g, '%23')
    .replace(/{/g, '%7B')
    .replace(/}/g, '%7D')
    .replace(/</g, '%3C')
    .replace(/>/g, '%3E')
    .replace(/atom/g, 'a')
    .replace(/bond/g, 'b')}") no-repeat center center`
}

self.onconnect = async function (e) {
  // initialize the rdkit
  if (!self.RDKit) {
    await initRDKitModule().then((res) => {
      self.RDKit = res
      self.RDKit.prefer_coordgen(true)
    }).catch((err) => {
      console.log(err)
    })
  }
  const port = e.ports[0]
  console.log('the ', i++, 'th time connected')
  port.onmessage = function (e) {
    props = JSON.parse(e.data)
    if (props) {
      // console.log('props',props)
      Promise.resolve(props)
        .then(data => cssBgSvg(initHighlight(renderMol(data))))
        .then((svg) => {
          if (props.css)
            port.postMessage(svg)

          else
            port.postMessage(svg)
        })
        .catch((err) => {
          console.log('promise in worker error:', err)
        })
      // console.log('out',out)
      // if (props.css){
        // port.postMessage(cssBgSvg(out));
      // }else{
        // port.postMessage(cssBgSvg(out));
      // }
      // console.log(cssBgSvg(out))
      // out=null;
    }
  }
  port.onmessageerror = function (e) {
    console.log('onmessageerror', e)
  }
}
self.onerror = function (e) {
  console.log('onerror', e)
}