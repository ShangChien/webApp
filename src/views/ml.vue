<script setup lang='ts'>
import mlBody from '@/components/ml/mlBody.vue'
</script>

<template>
  <div class="h-90vh relative m-0 p-0 flex ">
    <Suspense>
      <template #default>
        <mlBody />
      </template>
      <template #fallback>
        <span>Loading...</span>
      </template>
    </Suspense>
  </div>
</template>

<style>

</style>
