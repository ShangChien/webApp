export interface molData{
    smiles?:String
    qsmiles?:String
    atoms?:Number[]
    bonds?:Number[]
    addAtomIndices?:Boolean
    addBondIndices?:Boolean
    legend?:String
    width?:Number
    height?:Number
    highlightColor?:Number[]
    bondLineWidth?: Number
    highlightBondWidthMultiplier?:Number
    highlightRadius?:Number
    minFontSize?:Number
    explicitMethyl?:Boolean
}