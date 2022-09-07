export interface molData {
  id?: number;
  smiles?: string | any;
  qsmiles?: string;
  selected?: boolean;
  atoms?: {[key: number|string]: number[]};
  bonds?: {[key: number|string]: number[]};
  label?: (string|null)[];
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
