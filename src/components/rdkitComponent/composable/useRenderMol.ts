import type { molData } from "@/components/types";

export function useRenderMol(props: molData,rdkit:any) {
  var atomsIndex = Object.values(props.atoms??{}).reduce((pre,cur)=>pre.concat(cur),[])
  var bondsIndex = Object.values(props.bonds??{}).reduce((pre,cur)=>pre.concat(cur),[])
  //console.log(props.smiles,)
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
  mDetail["highlightRadius"] = props.highlightRadius ?? 0.3;
  mDetail["minFontSize"] = props.minFontSize ?? 10;
  mDetail = JSON.stringify(mDetail);
  let svg = mol.get_svg_with_highlights(mDetail);
  //console.log(svg)
  mol.delete();
  mol=null;
  qmol.delete();
  qmol=null;
  if (svg=='') {
    svg=`<svg><rect></rect></svg>`
    console.log('无效的smile式',props.smiles)
  }
  //console.log(svg)
  return svg;
}