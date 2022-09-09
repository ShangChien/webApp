import type { molData } from "@/components/types";

export function useRenderMol(props: molData,rdkit:any) {
  
  const concatByIndex = (list1:{[key: number|string]: number[]}|undefined) =>{ 
    let outList: any[] = []
    for(let i in list1){
      outList = outList.concat(list1[i])
    }
    return Array.from(new Set(outList)) 
  }
  var atomsIndex = concatByIndex(props.atoms)
  var bondsIndex = concatByIndex(props.bonds)
  console.log(props.smiles,)
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
  mol=null;
  qmol.delete();
  qmol=null;
  return svg;
}