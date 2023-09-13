<script setup lang='ts'>
import monaco from '@/components/monaco/fileSelector.vue'
</script>

<template>
  <div class="h-90vh relative m-0 p-0 flex ">
    <Suspense>
      <template #default>
        <monaco />
      </template>
      <template #fallback>
        <span>Loading...</span>
      </template>
    </Suspense>
  </div>
</template>

<style>

</style>
