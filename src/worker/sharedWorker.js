// import rdkit into the worker, and export the renderMol function
importScripts("/RDKit_minimal.js")
//initialize the rdkit
initRDKitModule().then(res=>{
	self.rdkit = res	
  self.rdkit.prefer_coordgen(true);
}).catch(err=>{
	console.log(err)
})

let i = 0 //counter for the number of times the worker is connected
let out = null //the output(svg image) of the renderMol function
function renderMol(props){
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
  out = mol.get_svg_with_highlights(mDetail);
  mDetail = null;
  mdetailsRaw = null;
  qmol.delete();
  qmol = null;
  mol.delete();
  mol = null;
  return out
}

self.onconnect = function(e) {
  var port = e.ports[0];
  console.log("the ",i++,"th time connected")
  port.onmessage = function(e) {
    if (e.data) {
      out = renderMol(JSON.parse(e.data))
      port.postMessage(out);
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






