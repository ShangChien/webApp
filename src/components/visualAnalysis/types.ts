export enum coordOption {
  cartesian = 'cartesian',
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
  view = 'view',
  spaceFlex = 'spaceFlex',
  facetRect = 'facetRect',
  repeatMatrix = 'repeatMatrix',
}

export enum scaleType {
  linear = 'linear',
  time = 'time',
  log = 'log',
}

export enum scaleAxisType {
  x = 'x',
  y = 'y',
  color = 'color',
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

export interface itemOption {
  type: fig
  data?: object[]
  encode: {
    x: string
    y: string | string[]
    color?: string
  }
  scale?: {
    x?: ScaleItem
    y?: ScaleItem
    color?: ScaleItem
  }
  style?: any
  shape?: any // customize
  tooltip?: any // customize
  coordinate?: {
    type: coordOption
    transform?: {
      value: { type: string }[] // 转置XY轴 [{ type: 'transpose' }]
    }
  }
}

export interface ViewType {
  type: view
  data: object[]
  children: itemOption[]
  coordinate?: {
    type: coordOption
    transform?: {
      value: { type: string }[] // 转置XY轴 [{ type: 'transpose' }]
    }
  }
  scale?: {
    x?: ScaleItem
    y?: ScaleItem
    color?: ScaleItem
  }
}
