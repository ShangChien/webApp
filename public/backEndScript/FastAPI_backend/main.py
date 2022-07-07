from fastapi import Depends, FastAPI
from fastapi.middleware.cors import CORSMiddleware
from dependencies import Smi_DummyAtom, Smi_DummyBond, handle_core, atom_enum
from pydantic import BaseModel
from typing import Union,List,Set
######处理ligand结构
# Smi_DummyAtom(smi='',atoms=[]):##[smi*,]根据smile和原子索引,输出带有*的smiles
# Smi_DummyBond(smi='',bonds=[]):##{'smi*':[dummy1,dummy2],}根据smile和bond索引,输出带有虚原子的smiles和裂键后虚原子的索引
######处理core结构
# handle_core(smi='',atoms=[],bonds=[]):##[smi,[atoms,],[[begin,end],]]处理core结构

class molecule(BaseModel):
    id: int
    smiles: str
    atoms: List[int]=[]
    bonds: List[int]=[]
    
class mol4enum(BaseModel):
    core: List[molecule]
    ligand: List[molecule]


app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/enum")
async def enum_molecule(mol4enum:mol4enum):
    dummy_atoms_ligand_smis=[]
    for ligand in mol4enum.ligand:
        dummy_atoms_ligand_smis.extend(Smi_DummyAtom(ligand.smiles,ligand.atoms))
    cores=[]
    for core in mol4enum.core:
        cores.append(handle_core(core.smiles,core.atoms,core.bonds))
    smiles_out=atom_enum(cores,dummy_atoms_ligand_smis)
    return {"message": smiles_out[0]}
    

