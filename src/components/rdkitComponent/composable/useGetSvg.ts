import { reactive } from "vue";
import type { molData } from "@/components/types";
import { useRenderMol } from '@/components/rdkitComponent/composable/useRenderMol'
const Color = (n:number) => 'hsla('+ Math.floor((n+8.6)*36) +',90%,70%,1)'
export function useGetSvg(props: molData,rdkit:any) {
	const svgItem:any= reactive({
		svg: {},
		rect: {},
		ellipse: [],
		path: {
			hightBonds: [],
			symble: [],
			bond: [],
		},
	});
	const svg=useRenderMol(props,rdkit)
	const parser = new DOMParser();
  let xmlDoc = parser.parseFromString(svg, "text/xml");
  //获取svgRoot内容
  let svgRoot = xmlDoc.getElementsByTagName("svg")[0];
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
    svgItem.svg[attrs] = svgRoot.getAttribute(attrs);
  }
  //获取svgRect内容
  let svgRect = xmlDoc.getElementsByTagName("rect")[0];
  for (let attrs of ["style", "width", "height", "x", "y"]) {
    svgItem.rect[attrs] = svgRect.getAttribute(attrs);
  }
  //获取svgEllipse内容（高亮无符号的原子）
  let svgEllipse: any = xmlDoc.getElementsByTagName("ellipse");
  for (var item of svgEllipse) {
    svgItem.ellipse.push({ ellipse: {} });
    for (let attrs of ["cx", "cy", "rx", "ry", "class"]) {
      svgItem.ellipse[svgItem.ellipse.length - 1].ellipse[attrs] =
        item.getAttribute(attrs);
    }
    svgItem.ellipse[svgItem.ellipse.length - 1].ellipse["style"] = item
      .getAttribute("style")
      .concat(";opacity: 0.4");
  }
  //获取svgPath内容
  let svgPath: any = xmlDoc.getElementsByTagName("path");
  for (let item of svgPath) {
    //字符匹配
    if (item.getAttribute("style") == null) {
      svgItem.path.symble.push({ path: {} });
      for (let attrs of ["class", "d", "fill"]) {
        svgItem.path.symble[svgItem.path.symble.length - 1].path[attrs] =
          item.getAttribute(attrs);
      }
      //化学键匹配
    } else if (
      item.getAttribute("style")?.split(";")[3] == "stroke-width:1.0px"
    ) {
      svgItem.path.bond.push({ path: {} });
      for (let attrs of ["class", "d", "style"]) {
        svgItem.path.bond[svgItem.path.bond.length - 1].path[attrs] =
          item.getAttribute(attrs);
      }
      //高亮化学键
    } else {
      svgItem.path.hightBonds.push({ path: {} });
      for (let attrs of ["class", "d", "style"]) {
        svgItem.path.hightBonds[svgItem.path.hightBonds.length - 1].path[
          attrs
        ] = item.getAttribute(attrs);
      }
      svgItem.path.hightBonds[svgItem.path.hightBonds.length - 1].path[
        "style"
      ] = item.getAttribute("style").concat(";opacity: 0.3");
    }
  }
  //如果props.atoms存在非空属性，则添加原子索引
  Object.keys(props?.atoms as {[key: number]: number[]}).forEach((key:string) => {
    props.atoms?.[key].forEach((n:number)=>{
      let item = svgItem.ellipse.findIndex((v:any)=>v.ellipse.class.split('-')[1] == n.toString())
      svgItem.ellipse[item].ellipse.style.replace('#9FACE6', Color(+key))
    })
  })
  

  console.log(svgItem)
	return svgItem
  
}