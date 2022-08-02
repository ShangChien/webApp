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
export interface mol4E {
  id: number
  smiles: string
  atoms?: (number|null)[]
  bonds?: (number|null)[]
  label: (string|null)[]
}
