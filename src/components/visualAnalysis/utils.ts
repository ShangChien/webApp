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
