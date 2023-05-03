from rdkit.Chem import rdMolEnumerator,rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit import Chem
from pycdxml import cdxml_slide_generator, cdxml_converter
from itertools import combinations,product
from joblib import Parallel, delayed
import copy

from pydantic import BaseModel
from typing import Union,Iterable,List,Set,Dict,Any

class enumSetting(BaseModel):
	list:List[int]=[]
	range:List[int]=[]
	connect2index:List[int]=[]
	keepSame2Index:List[int]=[]

class core(BaseModel):
	id: int
	smiles: str
	enumAtoms:Dict[int, enumSetting]={}
	enumBonds:Dict[int, enumSetting]={}

class molecule(BaseModel):
  id: int
  smiles: str
  atoms: Dict[int,List[int]]={}
  bonds: Dict[int,List[int]]={}
    
class enumData(BaseModel):
  core: core
  ligands: List[molecule]


def add_dummy2smi(mol:Any,i:int)->List[Union[str,int]]:
	dummy = Chem.Atom(0)  # 0 is the atomic number 
	rwmol = Chem.RWMol(mol)
	rwmol.AddBond(rwmol.AddAtom(dummy),i,Chem.rdchem.BondType.SINGLE)
	nums=rwmol.GetNumAtoms()
	smi=Chem.MolToSmiles(rwmol)
	return [smi,nums]

def split_bond2smi(mol:Any,i:int)->List[Union[str,int]]:
	ibond=mol.GetBondWithIdx(i)
	begin=ibond.GetBeginAtomIdx()
	end=ibond.GetEndAtomIdx()
	rwmol=Chem.RWMol(mol)
	rwmol.RemoveBond(begin,end)
	rwmol.ReplaceAtom(begin,Chem.Atom(0))
	rwmol.ReplaceAtom(end,Chem.Atom(0))
	smi=Chem.MolToSmiles(rwmol)
	nums=rwmol.GetNumAtoms()
	dummy_index=[
		atom.GetIdx()+1 
		for atom in Chem.MolFromSmiles(smi,sanitize=False).GetAtoms() 
			if atom.GetSymbol() == "*"
	]
	return [smi,nums,*dummy_index]

def pre_ligand(arr:List[molecule])->Dict[str,List[List]]:
	from rdkit import Chem
	colorDummyAtomsIndex={}
	colorDummyBondsIndex={}
	for ligand in arr:
		mol=Chem.MolFromSmiles(ligand.smiles)
		for color in ligand.atoms:
			if not (color in colorDummyAtomsIndex.keys()):
				colorDummyAtomsIndex[color]=[]
			colorDummyAtomsIndex[color].extend(
				list(map(lambda i:add_dummy2smi(mol,i),ligand.atoms[color]))
			)
		for color in ligand.bonds:
			if not (color in colorDummyBondsIndex.keys()):
				colorDummyBondsIndex[color]=[]
			colorDummyBondsIndex[color].extend(list(map(lambda i:split_bond2smi(mol,i),ligand.bonds[color])))
	return {'atoms':colorDummyAtomsIndex,'bonds':colorDummyBondsIndex}

def enum_core_site(core:core)->Iterable[Dict[int,List[int]]]:
	rateCombo={}
	for k in core.enumAtoms:
		arr=core.enumAtoms[k].list
		rate=sorted(core.enumAtoms[k].range)
		rate[-1]=rate[-1]+1
		rateCombo[k]=[]
		for n in range(*rate):
			rateCombo[k].extend([list(i) for i in combinations(arr, n)])
	combo=[[*i] for i in product(*rateCombo.values())]
	colorIndex=rateCombo.keys()
	combo=map((lambda x : {v:x[i] for i,v in enumerate(colorIndex)}), combo)
	return combo

def getAtomCombos(core:core,ligands:Dict[str,Dict[int,List]])->List[str]:
	#组合不同颜色的位点库
	mol=Chem.MolFromSmiles(core.smiles)
	nums=mol.GetNumAtoms()
	colorbase={}
	for i in core.enumAtoms.keys():
		colorbase[i]=[]
		for i1 in core.enumAtoms[i].connect2index:
			colorbase[i].extend(ligands['atoms'][i1])
	#组合带有虚原子的smi
	t_combos=[]
	combos=enum_core_site(core)
	for combo in combos:
		keys=combo.keys()
		keys_combo=[]
		for k in keys:
			keys_combo.append([])
			length=len(combo[k])
			for combo_i in product(colorbase[k],repeat=length):
				_item=[]#[[smi,totalAtoms,site]]
				for i,_v in enumerate(combo_i): 
					_convert=_v.copy()
					_convert.append(combo[k][i])
					_item.append(_convert)
				keys_combo[-1].append(_item)
		for combo_color in product(*keys_combo):
			i=[c for com in combo_color for c in com]
			smis=core.smiles
			site_n=nums
			link='m:'
			for item in i:
				smis=smis+'.'+item[0]
				link=link+str(site_n)+':'+str(item[-1])+','
				site_n=site_n+item[1]
			cxsmi=smis+' '+'|'+link+'|'
			if not (link=='m:'):
				t_combos.append(cxsmi) 
	return t_combos

def convert_smi(i:str)->str:
	mol=Chem.MolFromSmiles(i)
	smi=Chem.MolToSmiles(rdMolEnumerator.Enumerate(mol)[0])
	Chem.CanonSmiles(smi)
	return Chem.CanonSmiles(smi)

def parallel_convert_smis(smis_link:List[str])->List[str]:
	smis=Parallel(n_jobs=8)(delayed(convert_smi)(i) for i in smis_link)
	return smis

async def enum_atoms_smiles(enumData:enumData)->List[str]:
	ligands=pre_ligand(enumData.ligands)
	smi_link=getAtomCombos(core=enumData.core,ligands=ligands)
	smis_p= parallel_convert_smis(smi_link)
	smis=list(set(smis_p))
	return smis

def layout_cdxml(cdxmls:List[Any],cols:int,rows:int=-1,single_page:bool=True,scale:float=1.0)->List[Any]:
	##when single_page=false, must explicitly pass rows with value>0
	import math
	length = len(cdxmls)
	props = [[cdxml_slide_generator.TextProperty('index', i+1)] for i in range(length)]
	real_rows = math.ceil(length/cols) if single_page else rows
	#logic
	slides=[]
	if single_page:
		sg = cdxml_slide_generator.CDXMLSlideGenerator(
			style="ACS 1996", 
			number_of_properties=1, 
			columns=cols, rows=real_rows, 
			slide_width=cols*10*scale, 
			slide_height=real_rows*10*scale
		)
		slides.append(sg.generate_slide(cdxmls, props))
	else:
		if rows<0:
			print("Value Error: rows=", rows)
		else:
			per_page_nums=cols*rows
			total_pages=math.ceil(length/per_page_nums)
			def getSlide(pageIndex:int):
				##必须将初始化sg放到并行函数内部才能运行
				sg = cdxml_slide_generator.CDXMLSlideGenerator(
					style="ACS 1996", 
					number_of_properties=1, 
					columns=cols, rows=real_rows, 
					slide_width=cols*10*scale, 
					slide_height=real_rows*10*scale
				)
				start = pageIndex * per_page_nums
				end   = start + per_page_nums
				cdxml = cdxmls[start:end]
				prop = props[start:end]
				slide = sg.generate_slide(cdxml, prop)
				return slide
			#sub_cdxmls=[cdxmls[i*per_page_nums:(i+1)*per_page_nums] for i in range(total_pages)]
			#sub_props=[props[i*per_page_nums:(i+1)*per_page_nums] for i in range(total_pages)]
			#slides=[sg.generate_slide(sub_cdxmls[i], sub_props[i]) for i,_v in enumerate(sub_props)]
			slides=Parallel(n_jobs=10)(delayed(getSlide)(i) for i in range(total_pages))
	return slides