<script setup lang='ts'>
import { NTabPane, NTabs } from 'naive-ui'
import { useMlState } from '@/components/ml/mlState'

const { currentTab } = useMlState()
</script>

<template>
  <NTabs v-model:value="currentTab" type="segment">
    <NTabPane name="files" tab="Inspect file" />
    <NTabPane name="mlPredict" tab="Model predict" />
    <NTabPane name="mlFit" tab="Model fit" />
    <NTabPane name="visualAnalysis" tab="Visual analysis" />
  </NTabs>
</template>

<style>
</style>
