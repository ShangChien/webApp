<script setup lang="ts">
import { NButtonGroup, NButton, NGradientText,NPopover } from "naive-ui";
import { inject,computed,ref } from "vue";

const siteType:any= inject('siteType')
const options:any= Array.from(Array(10).keys()).map((i:number)=> {return {label:i,value:i}})         
const Color = computed(() => 'hsla('+(siteType.value+8.6)*36 +',90%,70%,1)')
const visiable=ref(false)
</script>
<template>
<div >
  <div class="i-fluent-emoji-flat-artist-palette text-3xl mb-1"></div>
  <div class="text-xl inline-block ml-1 ma-1" >: </div>
  <n-popover class="inline-block mb-1 cursor-pointer" trigger="manual" placement="top-start" :show="visiable" >
    <template #trigger>
      <div class="i-carbon-circle-filled text-3xl ma-1 rd-50% cursor-pointer outcircle 
                    hover:(transform-scale-110 transition-transform duration-210 ease-in-out)"
                    @click="visiable=!visiable" />
    </template>
    <div class="p-0 m--2 mb--1" @click="visiable=!visiable">
      <div v-for="option in options" 
        :key="option.value" 
        @click="siteType=option.value"
        :style="{color:'hsla('+(option.value+8.6)*36 +',90%,70%,1)'}"
        class="i-ic-sharp-circle mt-1 pt-1 text-3xl rd-0.5 
               hover:(i-carbon-circle-filled transform-scale-120 duration-210 ease-in-out)" />
    </div>
  </n-popover>
</div>
</template>
<style scoped>
.outcircle {
  color: v-bind('Color');
}
</style>
