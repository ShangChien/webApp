<script setup lang="ts">
import { ref, onMounted } from "vue";
import { NStep,NSteps,NCard,NButton,NButtonGroup,NIcon,NSpace,NDivider,NSwitch,NPageHeader,
NGrid,NGridItem,NStatistic,NAvatar,NDropdown,NThing } from "naive-ui"
import { MdArrowRoundBack,MdArrowRoundForward } from '@vicons/ionicons4'
import siteSelcet from "@/components/rdkitComponent/siteSelect.vue"

const currentRef=ref(0);
const showInputBoard=ref(true);
const options = [{ label: "Ketcher: 2.4.0" }, { label: "RDKit: 2022.3.2" }];
</script>

<template>
<div>
<n-page-header subtitle="分子结构枚举，子结构筛选 ( beta 1.0 )">
  <n-divider />
  <template #title>
    <a style="text-decoration: none; color: inherit">分子生成器</a>
  </template>
  <template #avatar>
    <n-avatar
      color=" white"
      src="/favicon.svg"
    />
  </template>
  <template #extra>
    <n-space>
      <n-switch v-model:value="showInputBoard" >
        <template #checked>
          已显示画板
        </template>
        <template #unchecked>
          已隐藏画板
        </template>
      </n-switch>
      <n-dropdown :options="options" placement="bottom-start">
        <n-button text>· · ·</n-button>
      </n-dropdown>
    </n-space>
  </template>
  <site-selcet v-show="showInputBoard" />
</n-page-header>

<n-card embedded>
  <n-steps v-model:current="currentRef">
    <n-step title="配体" description="标记配体位点和分类" />
    <n-step title="主体" description="标记主核位点和分类" />
    <n-step title="位点组合" description="位点组合" />
		<n-step title="结果展示" description="结果查看和筛选" />
  </n-steps>
	<n-space style="padding-top:10px" >
    <n-button @click="currentRef--" :disabled="currentRef in [0,1]" type="info">
      back
    </n-button>
		<n-button @click="currentRef++" :disabled="currentRef===4" type="info">
      {{currentRef === 0 ? 'initialze' : 'next'}}
    </n-button>
	</n-space>
</n-card>
<div v-if="currentRef === 1" style="padding-top:10px;">
  
</div>
<div v-else-if="currentRef === 2">
  主核
</div>
<div v-else-if="currentRef === 3">
  组合设置
</div>
<div v-else-if="currentRef === 4">
  结果页
</div>
<div v-else>
<n-thing>
  <n-grid :cols="5" >
    <n-grid-item>
      <n-statistic label="主核" value="125 个" />
    </n-grid-item>
    <n-grid-item>
      <n-statistic label="配体" value="22 位" />
    </n-grid-item>
    <n-grid-item>
      <n-statistic label="结构" value="36 次" />
    </n-grid-item>
    <n-grid-item>
      <n-statistic label="计算次数" value="83 个" />
    </n-grid-item>
  </n-grid>
</n-thing>
</div>
</div>
</template>