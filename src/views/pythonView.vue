<script setup lang='ts'>
import pythonComponent from '@/components/python/pythonComponent.vue'
</script>

<template>
  <div class="h-90vh relative m-0 p-0 flex ">
    <Suspense>
      <template #default>
        <python-component />
      </template>
      <template #fallback>
        <span>Loading...</span>
      </template>
    </Suspense>
  </div>
</template>

<style>

</style>
