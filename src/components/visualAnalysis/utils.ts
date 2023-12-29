export function readCsv(csvText: string): object[] {
  const lines = csvText.trim().split('\n')
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
