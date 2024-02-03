export function readCsv<T>(csvText: string): T[] {
  const lines = csvText.trim().split(/\r?\n/)
  const cols = lines[0].split(',').map(col => col.trim())
  const arr = []
  lines.slice(1).forEach((line, i) => {
    const items = line.split(',')
    if (items.length === cols.length) {
      const obj = items.reduce((obj, item, index) => {
        const trimmedItem = item.trim()
        const numberPattern = /^-?\d+(\.\d+)?$/
        obj[cols[index]] = numberPattern.test(trimmedItem) ? Number.parseFloat(trimmedItem) : trimmedItem
        return obj
      }, {})
      arr.push(obj)
    } else {
      console.error(`Line ${i + 1} does not have the same number of items as the header.`)
    }
  })
  return arr
}

export function file2Data<T>(file: { name: string; contents: string }): T[] {
  if (file === null || file === undefined || Object.keys(file).length === 0) {
    console.log(`${file} is empty!`)
    return []
  } else {
    let data: T[]
    if (file.name.endsWith('csv') || file.name.endsWith('CSV')) {
      data = readCsv<T>(file.contents)
    } else if (file.name.endsWith('json') || file.name.endsWith('JSON')) {
      data = JSON.parse(file.contents)
    }
    return data
  }
}

export function normalize<T>(data: T[]): T[] {
  if (data.length > 0) {
    const cols: string[] = Object.keys(data[0])
    cols.slice(1).forEach((col) => {
      const col_arr = data.map(item => item[col])
      const min = Math.min(...col_arr)
      const max = Math.max(...col_arr)
      data = data.map(item => ({
        ...item,
        [col]: (item[col] - min) / (max - min),
      }))
    })
    return data
  } else {
    console.log(`${data} is empty!`)
    return []
  }
}

export function toOptionsArray(enumValue: any): { label: string; value: any }[] {
  const options = []
  Object.keys(enumValue).forEach((key: string) => {
    options.push({
      label: key,
      value: enumValue[key],
    })
  })
  return options
}

export function tickInfo(min: number, max: number, step = 20): number[] {
  return Array.from({ length: Math.floor((max - min) / step) + 1 }, (_, index) => 250 + index * step)
}
