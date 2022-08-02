<script setup lang="ts">
import { defineComponent, ref } from 'vue'
import { repeat } from 'seemly'
import type { TreeOption } from 'naive-ui'
import { NTree } from 'naive-ui'

const props = defineProps();
const emit = defineEmits([]);

</script>

<template>
<n-tree 
	block-line
  cascade
  checkable
/>


</template>
<style></style>