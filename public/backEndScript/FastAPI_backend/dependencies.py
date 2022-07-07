from rdkit.Chem import rdMolEnumerator
from rdkit import Chem
from itertools import combinations

def Smi_DummyAtom(smi='',atoms=[]):##[smi,]
	smis=[]
	for i in atoms:
		mol=Chem.MolFromSmiles('**.'+smi+' |m:0:'+str(i+2)+'|')
		mol_=rdMolEnumerator.Enumerate(mol)[0]
		smis_=Chem.MolToSmiles(mol_)
		smis.append(smis_)
	return list(set(smis))

def Smi_DummyBond(smi='',bonds=[]):##{'smi':[bond],}
	smis={}
	mol=Chem.MolFromSmiles(smi)
	for i in bonds:
		begin=mol.GetBonds()[i].GetBeginAtomIdx()
		end=mol.GetBonds()[i].GetEndAtomIdx()
		mol1=Chem.RWMol(mol)
		mol1.RemoveBond(begin,end)
		mol1.ReplaceAtom(begin,Chem.Atom(0))
		mol1.ReplaceAtom(end,Chem.Atom(0))	
		mol_=mol1.GetMol()
		smis_dummy=Chem.MolToSmiles(mol_)
		#获取虚原子的索引
		dummy_atoms=[]
		for atom in Chem.MolFromSmiles(smis_dummy,False).GetAtoms():
			if atom.GetAtomicNum()==0:
				dummy_atoms.append(atom.GetIdx())
		smis[smis_dummy]=dummy_atoms
	return smis
	
def handle_core(smi='',atoms=[],bonds=[]):##处理core结构:[smi,[atoms,],[[begin,end],]]
	core_out=[]
	bonds_atomize=[]
	mol=Chem.MolFromSmiles(smi)
	smi_out=Chem.MolToSmiles(mol)
	for i in bonds:
		begin=mol.GetBonds()[i].GetBeginAtomIdx()
		end=mol.GetBonds()[i].GetEndAtomIdx()
		bonds_atomize.append([begin,end])
		bonds_atomize.append([end,begin])
	core_out=[smi_out,atoms,bonds_atomize]
	return core_out

##根据一个配体和主核进行atoms位点枚举
def atom_enum(core=[],dummy_atoms_ligand_smi=[],dummy_bonds_ligand_smi={}):##return:[[smi,],]
	smile_out_atoms=[]
	smile_out_bonds=[]
	smile4enum=[]
	for i in core:
		core_smi=i[0]
		core_atoms=i[1]
		core_bonds=i[2]
		##获取core_smi分子式的包含的原子个数
		core_atoms_num=Chem.MolFromSmiles(core_smi).GetNumAtoms()
		##获取core_smi分子式中被选中的虚原子个数
		dummy_atoms_num=len(core_atoms)
		for a_smi in dummy_atoms_ligand_smi:###????
			##获取dummy_atoms_ligand_smi分子式包含的原子个数
			ligand_atoms_num=Chem.MolFromSmiles(dummy_atoms_ligand_smi[0]).GetNumAtoms()
			smi_with_links=[]
			for substitute_num in range(1,dummy_atoms_num+1):
				smile_without_link=core_smi+('.'+a_smi)*substitute_num
				core_site_sets=combinations(core_atoms,substitute_num)
				dummy_site=[]
				for n in range(1,substitute_num+1):
					dummy_site.append(str(core_atoms_num+ligand_atoms_num*(n-1)))
				for core_site_set in core_site_sets:
					smi_link=' |m:'
					for n in range(1,substitute_num+1):
						smi_link=smi_link+str(dummy_site[n-1])+':'+str(core_site_set[n-1])+','
					smi_link=smi_link+'|'
					smi_with_links.append(smile_without_link+smi_link)
				smile4enum.extend(smi_with_links)
	for i in smile4enum:
		molbundle=rdMolEnumerator.Enumerate(Chem.MolFromSmiles(i))
		smile_out_atoms.append(Chem.MolToSmiles(molbundle.GetMol(0)))
	smile_out=[smile_out_atoms,smile_out_bonds]
	return smile_out