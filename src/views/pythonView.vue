<script setup lang='ts'>
import fileSelector from '@/components/monaco/fileSelector.vue'
</script>

<template>
  <div class="h-90vh relative m-0 p-0 flex ">
    <Suspense>
      <template #default>
        <fileSelector />
      </template>
      <template #fallback>
        <span>Loading...</span>
      </template>
    </Suspense>
  </div>
</template>

<style>

</style>
