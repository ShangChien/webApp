<script setup lang="ts">
import { NIcon, NButton,NDynamicTags, NGradientText,NPopselect } from "naive-ui";
import { Theater } from "@vicons/carbon";
import { h,inject } from "vue";

const molTags:any= inject('molTags')

</script>
<template>
	<n-dynamic-tags v-model:value="molTags"></n-dynamic-tags>
</template>