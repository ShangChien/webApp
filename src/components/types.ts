export interface molData {
  id?: number;
  name?: string;
  smiles?: string | any;
  qsmiles?: string;
  selected?: boolean;
  type?: 'core' | 'ligand'| 'mole';
  css?: boolean;
  atoms?: {[key: number]: number[]};
  bonds?: {[key: number]: number[]};
  labels?: (string|null)[];
  addAtomIndices?: boolean;
  addBondIndices?: boolean;
  legend?: string;
  width?: number;
  height?: number;
  highlightColor?: number[];
  bondLineWidth?: number;
  highlightBondWidthMultiplier?: number;
  highlightRadius?: number;
  minFontSize?: number;
  explicitMethyl?: boolean;
  editable?: boolean;
  extraData?: any;
}
export interface enumData {
	[key: number]:{
		list : number[],
		range : [number,number],
		keepSame2Index : number[],
		connect2index  : number[],
	}
}
