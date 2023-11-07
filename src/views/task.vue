<script setup lang="ts">
import {
  NButton,
  NButtonGroup,
  NIcon,
  NProgress,
  NSpace,
  NTable,
  NTbody,
  NTd,
  NTh,
  NThead,
  NTr,
} from 'naive-ui'
</script>

<template>
  <div>
    <NSpace vertical>
      <NTable :single-line="false" size="small">
        <NThead>
          <NTr>
            <NTh>序号</NTh>
            <NTh>名称</NTh>
            <NTh>类型</NTh>
            <NTh>状态</NTh>
            <NTh>操作</NTh>
          </NTr>
        </NThead>
        <NTbody>
          <NTr>
            <NTd>1</NTd>
            <NTd>220510-ETA-0028</NTd>
            <NTd>结构优化</NTd>
            <NTd>
              <NProgress
                type="line"
                :percentage="60"
                indicator-placement="inside"
                color="#52c234"
                processing
              />
            </NTd>
            <NTd>
              <NButtonGroup>
                <NButton size="small" round tertiary type="primary">
                  <template #icon>
                    <NIcon><div class="i-carbon-view" /></NIcon>
                  </template>
                  查看
                </NButton>
                <NButton size="small" round tertiary type="info">
                  <template #icon>
                    <NIcon><div class="i-carbon-edit" /></NIcon>
                  </template>
                  编辑
                </NButton>
                <NButton size="small" round tertiary type="error">
                  <template #icon>
                    <NIcon><div class="i-carbon-trash-can" /></NIcon>
                  </template>
                  删除
                </NButton>
              </NButtonGroup>
            </NTd>
          </NTr>
          <NTr>
            <NTd>...</NTd>
            <NTd>...</NTd>
            <NTd>...</NTd>
            <NTd>...</NTd>
            <NTd>...</NTd>
          </NTr>
          <NTr>
            <NTd>...</NTd>
            <NTd>...</NTd>
            <NTd>...</NTd>
            <NTd>...</NTd>
            <NTd>...</NTd>
          </NTr>
        </NTbody>
      </NTable>
    </NSpace>
  </div>
</template>
