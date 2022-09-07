<script setup lang="ts">
import { ref, onMounted, watch, reactive, onUpdated,inject } from "vue";
import type { molData } from "@/components/types";
const rdkit:any = inject("rdkit");
const props = defineProps<molData>();
const emit = defineEmits(["update-mol"]);
const highlightMap= reactive({
  id: 0,
  smiles: props.smiles,
  atoms: {},
  bonds: {},
  label: [],
});
const svgitem: any = reactive({
  svg: {},
  rect: {},
  ellipse: [],
  path: {
    hightBonds: [],
    symble: [],
    bond: [],
  },
});
const parser = new DOMParser();
const xmlDoc = ref<Document>();
const siteType:any= inject('siteType')
const Color = (n:any) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,70%,1)'
function getSvgData(svg: string) {
  svgitem.svg = {};
  svgitem.rect = {};
  svgitem.ellipse = [];
  svgitem.path = {
    hightBonds: [],
    symble: [],
    bond: [],
  };
  xmlDoc.value = parser.parseFromString(svg, "text/xml");
  //获取svgRoot内容
  let svgRoot = xmlDoc.value.getElementsByTagName("svg")[0];
  for (var attrs of [
    "xmlns",
    "xmlns:rdkit",
    "xmlns:xlink",
    "version",
    "baseProfile",
    "xml:space",
    "width",
    "height",
    "viewBox",
  ]) {
    svgitem.svg[attrs] = svgRoot.getAttribute(attrs);
  }
  //获取svgRect内容
  let svgRect = xmlDoc.value.getElementsByTagName("rect")[0];
  for (let attrs of ["style", "width", "height", "x", "y"]) {
    svgitem.rect[attrs] = svgRect.getAttribute(attrs);
  }
  //获取svgEllipse内容（高亮无符号的原子）
  let svgEllipse: any = xmlDoc.value.getElementsByTagName("ellipse");
  for (var item of svgEllipse) {
    svgitem.ellipse.push({ ellipse: {} });
    for (let attrs of ["cx", "cy", "rx", "ry", "class"]) {
      svgitem.ellipse[svgitem.ellipse.length - 1].ellipse[attrs] =
        item.getAttribute(attrs);
    }
    svgitem.ellipse[svgitem.ellipse.length - 1].ellipse["style"] = item
      .getAttribute("style")
      .concat(";opacity: 0.6");
  }
  //获取svgPath内容
  let svgPath: any = xmlDoc.value.getElementsByTagName("path");
  for (let item of svgPath) {
    //字符匹配
    if (item.getAttribute("style") == null) {
      svgitem.path.symble.push({ path: {} });
      for (let attrs of ["class", "d", "fill"]) {
        svgitem.path.symble[svgitem.path.symble.length - 1].path[attrs] =
          item.getAttribute(attrs);
      }
      //化学键匹配
    } else if (
      item.getAttribute("style")?.split(";")[3] == "stroke-width:1.0px"
    ) {
      svgitem.path.bond.push({ path: {} });
      for (let attrs of ["class", "d", "style"]) {
        svgitem.path.bond[svgitem.path.bond.length - 1].path[attrs] =
          item.getAttribute(attrs);
      }
      //高亮化学键
    } else {
      svgitem.path.hightBonds.push({ path: {} });
      for (let attrs of ["class", "d", "style"]) {
        svgitem.path.hightBonds[svgitem.path.hightBonds.length - 1].path[
          attrs
        ] = item.getAttribute(attrs);
      }
      svgitem.path.hightBonds[svgitem.path.hightBonds.length - 1].path[
        "style"
      ] = item.getAttribute("style").concat(";opacity: 0.2");
    }
  }
}

function renderMol(props: molData) {
  const concatIndex = (list1:{[key: number|string]: number[]}|undefined) =>{ 
    let outList: any[] = []
    for(let i in list1){
      outList = outList.concat(list1[i])
    }
    return Array.from(new Set(outList)) 
  }
  var atomsIndex = concatIndex(props.atoms)
  var bondsIndex = concatIndex(props.bonds)
  let mol = rdkit.get_mol(props.smiles ?? "");
  let qmol = rdkit.get_qmol(props.qsmiles ?? "");
  //let mDetail = JSON.parse(mol.get_substruct_match(qmol))
  let mdetailsRaw = mol.get_substruct_matches(qmol);
  let mDetail = mdetailsRaw.length > 2 ? JSON.parse(mdetailsRaw) : [];
  mDetail = mDetail.reduce(
    (acc: any, { atoms, bonds }: any) => ({
      atoms: [...acc.atoms, ...atoms],
      bonds: [...acc.bonds, ...bonds],
    }),
    { bonds: [], atoms: [] }
  );
  mDetail["atoms"] = mDetail.atoms?.concat(atomsIndex ?? []) ?? atomsIndex;
  mDetail["bonds"] = mDetail.bonds?.concat(bondsIndex ?? []) ?? bondsIndex;
  mDetail["addAtomIndices"] = props.addAtomIndices;
  mDetail["addBondIndices"] = props.addBondIndices;
  mDetail["width"] = props.width ?? 200;
  mDetail["height"] = props.height ?? 200;
  mDetail["highlightColour"] = props.highlightColor ?? [0.624, 0.675, 0.902];
  mDetail["bondLineWidth"] = props.bondLineWidth ?? 1;
  mDetail["highlightBondWidthMultiplier"] =
    props.highlightBondWidthMultiplier ?? 20;
  mDetail["highlightRadius"] = props.highlightRadius ?? 0.25;
  mDetail["minFontSize"] = props.minFontSize ?? 10;
  //console.log(mDetail)
  mDetail = JSON.stringify(mDetail);
  let svg = mol.get_svg_with_highlights(mDetail);
  mol.delete();
  qmol.delete();
  getSvgData(svg);
  //console.log(svg)
  //console.log(bondMap.value)
  //console.log(svgitem);
  //rdkitdiv.value.innerHTML=svg
}
//目的,遍历highlightMap对象atoms属性的所有子属性，如果包含atomIndex，则删除
function preHandleIndex(obj:object|undefined,atomIndex:number){
    for (let key in obj) {
      if (obj[key].includes(atomIndex)) {
        obj[key].splice(
          obj[key].indexOf(atomIndex),
          1
        );
      }
    }
}
function domClick($event: any) {
  let itemList = $event.target.getAttribute("class").split(" ")[0].split("-");
  if (itemList[0] == "atom") {
    //如果点击atom
    //改颜色
    $event.target.style.stroke = Color(siteType.value);
    $event.target.style.fill = Color(siteType.value);
    $event.target.style.opacity = 0.8;
    var atomIndex = itemList[1] / 1
    preHandleIndex(highlightMap.atoms,atomIndex)
    //添加index到数组
    highlightMap.atoms[siteType.value] = highlightMap.atoms[siteType.value] ? 
                                         highlightMap.atoms[siteType.value]:[]
    highlightMap?.atoms[siteType.value].push(itemList[1] / 1);
    //highlightMap.atoms[siteType.value] = Array.from(new Set(highlightMap.atoms[siteType.value])).sort();
    emit("update-mol", highlightMap);
    //console.log('ssddd',highlightMap.atoms)
  } else if (itemList[0] == "bond") {
    //如果点击bond
    $event.target.style.stroke = Color(siteType.value);
    $event.target.style.fill = Color(siteType.value);
    $event.target.style.opacity = 0.4;
    var bondIndex = itemList[1] / 1
    preHandleIndex(highlightMap.atoms,bondIndex)
    //添加index到数组
    highlightMap.bonds[siteType.value] = highlightMap.bonds[siteType.value] ? 
                                         highlightMap.bonds[siteType.value]:[]
    highlightMap?.bonds[siteType.value].push(itemList[1] / 1);
    //highlightMap.bonds[siteType.value] = Array.from(new Set(highlightMap.bonds[siteType.value])).sort();
    emit("update-mol", highlightMap);
    //console.log(highlightMap.bonds)
  } else {
    console.log(itemList[0], "error");
  }
}
function domDblClick($event: any) {
  let itemList = $event.target.getAttribute("class").split(" ")[0].split("-");
  if (itemList[0] == "atom") {
    //如果点击atom
    //改颜色
    $event.target.style.stroke = "#9FACE6";
    $event.target.style.fill = "#9FACE6";
    $event.target.style.opacity = 0.6;
    //删除index
    preHandleIndex(highlightMap.atoms,itemList[1] / 1)
    emit("update-mol", highlightMap);
    //console.log(highlightMap.highlightAtoms)
  } else if (itemList[0] == "bond") {
    //如果点击bond
    $event.target.style.stroke = "#9FACE6";
    $event.target.style.fill = "#9FACE6";
    $event.target.style.opacity = 0.2;
    //删除index
    preHandleIndex(highlightMap.bonds,itemList[1] / 1)
    emit("update-mol", highlightMap);
    //console.log(highlightMap.highlightAtoms)
  } else {
    console.log(itemList[0], "error");
  }
}
function clearAll() {
  renderMol(props);
  highlightMap.atoms = {};
  highlightMap.bonds = {};
  emit("update-mol", highlightMap);
}

//根据highlightMap初始化高亮svg中的atoms和bonds
const svg_id = ref();
// function initHighlightSvg() {
//   //初始化处理高亮atoms
//   for (let item of svg_id.value.getElementsByTagName("ellipse")) {
//     //在atoms中添加item.class属性具有highlightMap.atoms的元素
//     if (
//       highlightMap.atoms?.includes(item.getAttribute("class").split("-")[1] / 1)
//     ) {
//       item.style.stroke = Color(siteType.value);
//       item.style.fill = Color(siteType.value);
//       item.style.opacity = 0.8;
//     }
//   }
//   //初始化处理高亮bonds
//   for (let item of svg_id.value.getElementsByTagName("path")) {
//     if (
//       item.style.strokeWidth === "13.7px" &&
//       highlightMap.bonds?.includes(
//         item.getAttribute("class").split(" ")[0].split("-")[1] / 1
//       )
//     ) {
//       item.style.stroke = Color(siteType.value);
//       item.style.fill = Color(siteType.value);
//       item.style.opacity = 0.4;
//     }
//   }
// }

onMounted(async () => {
  rdkit.prefer_coordgen(true);
  renderMol(props);
  //console.log("onmounted",props)
  //initHighlightSvg()
});
// onUpdated(() => {
//   initHighlightSvg();
// });

watch(props, (newVal) => {
  renderMol(newVal);
  highlightMap.smiles = newVal.smiles;
  highlightMap.atoms = {};
  highlightMap.bonds = {};
  console.log("watch", highlightMap.atoms);
});
</script>

<template class="svg">
  <svg
    v-bind="svgitem.svg"
    ref="svg_id"
    class="svgstyle"
    @click.right.prevent="clearAll"
    style="width: 100%; height: 100%"
  >
    <rect v-bind="svgitem.rect" />
    <path
      v-for="item in svgitem.path.symble"
      v-bind="item.path"
    />
    <path
      v-for="item in svgitem.path.bond"
      v-bind="item.path"
    />
    <path
      v-for="item in svgitem.path.hightBonds"
      v-bind="item.path"
      @click="domClick($event)"
      @dblclick="domDblClick($event)"
    />
    <ellipse
      v-for="item in svgitem.ellipse"
      v-bind="item.ellipse"
      @click="domClick($event)"
      @dblclick="domDblClick($event)"
    />
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
