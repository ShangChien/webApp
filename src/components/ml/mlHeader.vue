<script setup lang='ts'>
import { inject } from 'vue'
import { NTabPane, NTabs } from 'naive-ui'

const headerTabName = inject('headerTabName', 'files')
</script>

<template>
  <NTabs v-model:value="headerTabName" type="segment">
    <NTabPane name="files" tab="inspect file" />
    <NTabPane name="ml" tab="ml flow" />
  </NTabs>
</template>

<style>
</style>
