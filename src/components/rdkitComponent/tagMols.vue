<script setup lang="ts">
import { NDynamicTags } from "naive-ui";
import { inject } from "vue";

const molTags:any= inject('molTags')

</script>
<template>
	<n-dynamic-tags type="info" v-model:value="molTags"></n-dynamic-tags>
</template>