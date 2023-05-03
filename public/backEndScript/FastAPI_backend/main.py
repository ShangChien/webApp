from re import L
import uvicorn
from joblib import Parallel,delayed
from rdkit import Chem
from fastapi import Depends, FastAPI,Body
from fastapi.middleware.cors import CORSMiddleware
from dependencies import enumData,enum_atoms_smiles,layout_cdxml
from pycdxml import cdxml_slide_generator, cdxml_converter

from pydantic import BaseModel
from typing import Union,Iterable,List,Set,Dict,Any

class layout(BaseModel):
  cols:int=10
  rows:int=20
  scale:float=1.0
  singlePage:bool=False

class smi2cdx(BaseModel):
  smiles:List[str]=[]
  layout:layout

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
async def convert2cdxml(data:smi2cdx=Body()):
  def smi2cdxml(smi):
    mol = Chem.MolFromSmiles(smi)
    return cdxml_converter.mol_to_document(mol).to_cdxml()
  cdxmls=Parallel(n_jobs=10)(delayed(smi2cdxml)(smi) for smi in data.smiles)
  slides=layout_cdxml(cdxmls = cdxmls,
                        cols = data.layout.cols,
                        rows = data.layout.rows,
                       scale = data.layout.scale,
                 single_page = data.layout.singlePage)
  return {'data':slides}
    
if __name__ == "__main__":
  uvicorn.run(app="main:app", 
              host="0.0.0.0", 
              port=8000, 
              reload=True, 
              debug=True)
