<script setup lang='ts'>
import { NTabPane, NTabs } from 'naive-ui'
import { useDataCheckStore } from '@/components/dataCheck/dataCheckStore'

const MLstore = useDataCheckStore()
</script>

<template>
  <NTabs v-model:value="MLstore.currentTab" type="segment">
    <NTabPane name="files" tab="文本内容" />
    <NTabPane name="spectrum" tab="光谱比对" />
    <NTabPane name="table" tab="表格比对" />
  </NTabs>
</template>

<style>
</style>
