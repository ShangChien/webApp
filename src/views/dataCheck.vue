<script setup lang='ts'>
import dataCheckBody from '@/components/dataCheck/dataCheckBody.vue'
</script>

<template>
  <div class="h-90vh relative m-0 p-0 flex box-border flex-auto">
    <Suspense>
      <template #default>
        <dataCheckBody />
      </template>
      <template #fallback>
        <span>Loading...</span>
      </template>
    </Suspense>
  </div>
</template>

<style>

</style>
