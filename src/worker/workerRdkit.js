// import rdkit into the worker, and export the renderMol function
importScripts("./RDKit_minimal.js")

function renderMol(props){
  self.rdkit.prefer_coordgen(true);
  let mol = self.rdkit.get_mol(props.smiles ?? "");
  let qmol = self.rdkit.get_qmol(props.qsmiles ?? "");
  let mdetailsRaw = mol.get_substruct_matches(qmol);
  let mDetail = mdetailsRaw.length > 2 ? JSON.parse(mdetailsRaw) : [];
  mDetail = mDetail.reduce(
    (acc, { atoms, bonds }) => ({
      atoms: [...acc.atoms, ...atoms],
      bonds: [...acc.bonds, ...bonds],
    }),
    { bonds: [], atoms: [] }
  );
  mDetail["atoms"] = mDetail.atoms?.concat(props.atoms ?? []) ?? props.atoms;
  mDetail["bonds"] = mDetail.bonds?.concat(props.bonds ?? []) ?? props.bonds;
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
  let out = mol.get_svg_with_highlights(mDetail);
  qmol.delete();
  mol.delete();
  return out
}

onmessage = function (e) {
  initRDKitModule()
	.then(res=>{
		self.rdkit = res
    let out = renderMol(e.data)
    //console.log(out)
    sendMessage(out)
    self.rdkit=null
    self.close()
	})
	.catch(err=>{
		console.log(err)
    self.rdkit=null
    self.close()
	})
}






