export function figOptions(figType: string) {
  const encode = ['x', 'y', 'color', 'shape', 'size', 'opacity', 'series']

  const scaleOptions = {
    type: ['linear', 'log', 'time', 'ordina', 'point', 'quantile', 'quantize', 'threshold'],
    independent: [true, false],
    nice: [true, false],
    domain: [10, 100],
    range: [0, 1],
    key: '',
  }
  const scale = [
    { x: {} },
    { y: {} },
    { color: { type: [], range: [] } },
  ]
  // 分装G2的图形配置项，根据不同的mark类型，设置不同的配置项
  const options = {
    line: {
      encode,
    },
    bar: {},
    pie: {},
    // ...其他图表类型
  }
  return { type: figType, ...options[figType] }
}
