// uno.config.ts
import {
  defineConfig, presetIcons, presetUno,transformerVariantGroup
} from 'unocss'

export default defineConfig({
  presets: [
    presetUno(),
    presetIcons(),
  ],
  transformers: [
    transformerVariantGroup(),
  ],
  rules: [
    // your custom rules
    [
      /^hsla-(\d+)-(\d+)-(\d+)-(\d+)$/, 
      ([_,a,b,c,d]:any) => ({'background-color': `hsla(${a}, ${b}%, ${c}%, ${d}%)`})
    ],
    [
      /^font-([a-zA-Z]+)-([a-zA-Z]+)$/,
      ([_,a,b]:any)=>({'font-family': `${a}-${b}` })
    ]
  ],
})