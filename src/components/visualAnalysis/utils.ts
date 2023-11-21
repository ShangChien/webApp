export function _readCsv(csvText: string): object[] {
  const lines = csvText.split('\n')
  const arr = []

  for (const line of lines) {
    const obj = {}
    const items = line.split(',')
    for (const item of items) {
      obj[item.trim()] = item.trim()
    }
    arr.push(obj)
  }

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
