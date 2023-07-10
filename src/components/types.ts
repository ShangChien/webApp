import type { Component, InjectionKey, Ref, VNode } from 'vue'

export interface molData {
  id?: number
  name?: string
  smiles?: string | any
  qsmiles?: string
  selected?: boolean
  type?: 'core' | 'ligand' | 'mole'
  css?: boolean
  atoms?: { [key: number]: number[] }
  bonds?: { [key: number]: number[] }
  labels?: (string | null)[]
  addAtomIndices?: boolean
  addBondIndices?: boolean
  legend?: string
  width?: number
  height?: number
  highlightColor?: number[]
  bondLineWidth?: number
  highlightBondWidthMultiplier?: number
  highlightRadius?: number
  minFontSize?: number
  explicitMethyl?: boolean
  editable?: boolean
  extraData?: any
}
export interface enumData {
  [key: number]: {
    list: number[]
    range: [number, number]
    keepSame2Index: number[]
    connect2index: number[]
  }
}
interface electronCross {
  energy: number[]
  oscillator: number[]
}

export interface record {
  id: number
  timestamp: number
  length: number
  conditions: condition[]
}
export interface recordFull extends record {
  res: {
    id: number
    smiles: string
    name_calc: string
    name_mat: string
  }[]
}

export interface pgData {
  id: number
  smiles: string
}
export interface pgDataItem {
  id?: number
  timestamp?: number
  name_calc?: string
  name_mat?: string
  smiles?: string
  functional?: string
  basisset?: string
  charge?: number
  spinmultiplicity?: number
  numberbasisfunc?: number
  elections?: number[]
  ispcm?: boolean
  isspin?: boolean
  stationarytype?: string
  extrainfo?: string
  homo?: number
  lumo?: number
  eg?: number
  polar?: number
  energy?: number
  sfp?: any
  bfp?: any
  dirpath?: {
    ground?: string
    udft?: string
    tdft?: string
    polar?: string
    uv?: string
  }
  corrections?: {
    zeroPoint?: number
    energy?: number
    enthalpy?: number
    gibbsFreeEnergy?: number
  }
  structure?: {
    ground?: string
    udft?: string
    tdft?: string
  }
  tdft?: {
    s0_s1?: electronCross
    s0_t1?: electronCross
  }
  udft?: {
    energy?: number
    correction?: {
      zeroPoint?: number
      energy?: number
      enthalpy?: number
      gibbsFreeEnergy?: number
    }
  }
  spectrum?: {
    uv?: electronCross
    uv_low?: electronCross
  }
}

enum logicType {
  And = 'and',
  Not = 'not',
  Or = 'or',
}

export interface option {
  label?: string
  editable?: boolean
  value?: string
  key?: string
  icon?: VNode | Component
  disabled?: boolean
  items?: option[]
  children?: option[]
}

export interface condition extends Omit<option, 'value'> {
  id: number
  logic: logicType
  label?: string
  type?: string
  value?: string | number | string[] | number[] | Ref<number>[] | Ref<string>[] | Ref | any
  valueFormat?: string | number | string[] | number[] | Ref<number>[] | Ref<string>[] | Ref | any
  showDetail?: Ref<boolean>[] | Ref
  logic_icon?: Component
  label_icon?: VNode | Component
  component?: Component
}
export const keyStateKetcher: InjectionKey<{
  id: Ref<number>
  showModal: Ref<boolean>
  smiles: Ref<string>
}> = Symbol('stateKetcher')
