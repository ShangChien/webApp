import uvicorn
from joblib import Parallel,delayed
from rdkit import Chem
from fastapi import Depends, FastAPI,Body
from fastapi.middleware.cors import CORSMiddleware
from dependencies import enumData,enum_atoms_smiles,layout_cdxml
from pycdxml import cdxml_slide_generator, cdxml_converter
from typing import Union,Iterable,List,Set,Dict,Any

app = FastAPI()
app.add_middleware(
  CORSMiddleware,
  allow_origins=["*"],
  allow_credentials=True,
  allow_methods=["*"],
  allow_headers=["*"],
)

@app.post("/enum")
async def enum_molecule(enumData:enumData=Body()):
	smis=await enum_atoms_smiles(enumData)
	return {
    'message':'hello,just get it! ',
    'data':smis
  }
@app.post("/cdxml")
async def enum_molecule(smis:List[str]=Body()):
  def smi2cdxml(smi):
    mol = Chem.MolFromSmiles(smi)
    return cdxml_converter.mol_to_document(mol).to_cdxml()
  cdxmls=Parallel(n_jobs=6)(delayed(smi2cdxml)(smi) for smi in smis)
  slides=layout_cdxml(cdxmls=cdxmls,cols=6,single_page=False,rows=30)
  return {'data':slides}
    
if __name__ == "__main__":
  uvicorn.run(app="main:app", 
              host="0.0.0.0", 
              port=8000, 
              reload=True, 
              debug=True)
