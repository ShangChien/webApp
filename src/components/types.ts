export interface molData{
    id?: number
    smiles?:string|any
    qsmiles?:string
    atoms?:number[]
    bonds?:number[]
    addAtomIndices?:boolean
    addBondIndices?:boolean
    legend?:string
    width?:number
    height?:number
    highlightColor?:number[]
    bondLineWidth?: number
    highlightBondWidthMultiplier?:number
    highlightRadius?:number
    minFontSize?:number
    explicitMethyl?:boolean
    editable?:boolean
}