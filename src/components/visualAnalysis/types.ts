export enum coordOption {
  default = null,
  theta = 'theta',
  polar = 'polar',
  radial = 'radial',
}

export enum fig {
  interval = 'interval',
  bin = 'bin',
  line = 'line',
  point = 'point',
  area = 'area',
  pie = 'pie',
  cell = 'cell',
  rose = 'rose',
  radial = 'radial',
  radar = 'radar',
  box = 'box',
  boxplot = 'boxplot',
  density = 'density',
  heatmap = 'heatmap',
}

export enum view {
  default = 'view',
  spaceFlex = 'spaceFlex',
  facetRect = 'facetRect',
  repeatMatrix = 'repeatMatrix',
}

export interface Coord {
  type: coordOption
  transform?: { type: string }[] // 转置XY轴 [{ type: 'transpose' }]
}

export interface ScaleItem {
  type?: string
  range?: [number, number] | string[]
  domain?: [number, number]
  independent?: boolean
  key?: string
  nice?: boolean
  padding?: number
}

export interface Scale {
  x?: ScaleItem
  y?: ScaleItem
  color?: ScaleItem
}

export interface Option {
  type: fig
  data?: object[]
  encode: {
    x: string
    y: string | string[]
    color?: string
  }
  scale?: Scale
  style?: any
  shape?: any // customize
  tooltip?: any // customize
  coordinate?: Coord
}

export interface ViewType {
  type: view
  data: object[]
  children: Option[]
  coordinate?: Coord
  scale?: Scale
}
