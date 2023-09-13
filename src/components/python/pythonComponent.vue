<script setup lang='ts'>
import { onMounted, ref } from 'vue'
import { useElementSize, watchDebounced } from '@vueuse/core'
import { Pane, Splitpanes } from 'splitpanes'
import MonacoEditor from 'monaco-editor-vue3'
import { NButton, NScrollbar } from 'naive-ui'

const python = await import('@/worker/pyodide' /* @vite-ignore */)
// python setting
const py = ref<any>()
const codeStr = ref<string>(`
arr1=[1,2,3,4,5,6,7]
arr2=[x*x for x in arr1]
arr2
`)
const res = ref<string>(`截至 2022.10.9，Github发布上一些开源中文字体，均为由开源日文字体转化。现将这些字体列出如下。（先按授权方式，再按发行时间排序）\n
对于 SIL Open Font License 1.1 授权字体的使用许可与限制说明如下。 （适用于「小赖字体」「悠哉字体」「漫黑」「文楷」「975 系列字体」）

允许做的事：
这款字体无论是个人还是企业都可以自由商用，无需付费，也无需知会或者标明原作者。
这款字体可以自由传播、分享，或者将字体安装于系统、软件或APP中也是允许的，可以与任何软件捆绑再分发或一并销售。
这款字体可以自由修改、改造，制作衍生字体。修改或改造后的字体也必须同样以 SIL OFL 公开。
注意事项：
在制作衍生字体时，字体名称不可使用原有字体的「保留名称」 （如思源黑体衍生的字体不能以「思源」「Source」等名称命名，具体见各字体的授权说明） 。
根据 SIL Open Font License 1.1 的规定，禁止单独出售字体文件(OTF/TTF文件)的行为。
简言之：
可以用、可以嵌、可以改、不能单独卖。
`)
async function init() {
  if (!py.value) {
    console.log('init python env...')
    py.value = await python.loadPyodide()
    await py.value.loadPackage('micropip', { checkIntegrity: false })
    const micropip = py.value.pyimport('micropip')
    await micropip.install('numpy')
    py.value.runPython(`import sys
import numpy as np
print('python.version:', sys.version)
print('numpy.version:', np.__version__)`)
  } else {
    py.value = null
    console.log('python env destroyed')
  }
}
function run(code: string): void {
  res.value = py.value.runPython(code)
}
// monaco editor setting
const editorDom = ref(null)
const renderKey = ref<number>(0)
const { width: editorW, height: editorH } = useElementSize(editorDom)
const options = {
  colorDecorators: true,
  lineHeight: 2,
  tabSize: 2,
  theme: 'vs',
  language: 'python',
}
const operation = ['open', 'edit', 'view', 'about']
onMounted(() => {
  init()
})
watchDebounced(
  [editorH, editorW],
  () => {
    renderKey.value += 1
    console.log('key', editorW.value, editorH.value, renderKey.value)
  },
  { debounce: 1000, maxWait: 5000 },
)
</script>

<template>
  <Splitpanes horizontal class="splitpanes flex-auto relative" style="background-color: #ffffff">
    <Pane size="70" min-size="20">
      <div ref="editorDom" class="h-full flex flex-nowrap flex-col bg-slate-1 rd-1 min-w-280px">
        <div class="flex-none flex flex-nowrap justify-between ">
          <div class="flex flex-nowrap flex-none justify-start items-center m-0 p-0">
            <div v-for="(item, index) in operation" :key="index" class="rd-1 m-1 p-1 bg-slate-2">
              {{ item }}
            </div>
          </div>
          <div class="flex flex-nowrap flex-none justify-end items-center ">
            <div
              class="i-simple-icons-python text-1.5xl m-1 p-1"
              :class="[py ? 'bgColor' : '']"
              @click="init"
            />
            <NButton class="m-1" size="small" :disabled="!py" type="info" @click="run(codeStr)">
              run
            </NButton>
          </div>
        </div>
        <div :key="renderKey" v-eslint-disable class="pt-0 m-1 mt-0 flex-auto">
          <MonacoEditor
            v-model:value="res"
            :options="options"
            :height="editorH - 43"
            :weight="editorW"
          />
        </div>
      </div>
    </Pane>
    <Pane size="30" min-size="10">
      <div class="h-full mt-0 bg-slate-1 rd-1 flex flex-nowrap flex-col">
        <div class="bg-slate-2 rd-1 p-1 text-xl flex-none">
          Output:
        </div>
        <NScrollbar class="flex-auto text-2xl">
          <div>
            {{ res }}
          </div>
        </NScrollbar>
      </div>
    </Pane>
  </Splitpanes>
</template>

<style>
@import "splitpanes/dist/splitpanes.css";
.bgColor {
background-image:linear-gradient(314deg, #00c3ff, #ffff1c)
}

.splitpanes {background-color: #f8f8f8;}

.splitpanes__splitter {
  background-color:#6bd8fca3;
  position: relative;
  border-radius: 0.2rem;
  margin: 0.2rem;
}
.splitpanes__splitter:before {
  content: '';
  position: absolute;
  left: 0;
  top: 0;
  border-radius: 0.4rem;
  transition: opacity 0.4s;
  background-color: #6bd8fca3;
  opacity: 0;
  z-index: 1;
}
.splitpanes__splitter:hover:before {opacity: 1;}
.splitpanes--horizontal > .splitpanes__splitter:before { top:-2px;bottom:-2px;width: 100%;}
</style>
