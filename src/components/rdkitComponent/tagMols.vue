<script setup lang="ts">
import { storeToRefs } from 'pinia'
import { NDynamicTags } from "naive-ui";
import { inject } from "vue";
import { useEnumStore } from '@/stores/enumStore'
const enumStore = useEnumStore()
const { labels:molTags } =enumStore

</script>
<template>
	<n-dynamic-tags type="info" :value="molTags"></n-dynamic-tags>
</template>