export interface molData {
  id?: number;
  smiles?: string | any;
  qsmiles?: string;
  atoms?: number[] | any;
  bonds?: number[] | any;
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
export interface class4AtomsAndBonds {
  [key: number|string]: number[];
}
export interface mol4E {
  id: number
  smiles: string
  atoms?: class4AtomsAndBonds
  bonds?: class4AtomsAndBonds
  label: (string|null)[]
}
