
<script setup lang="ts">
import HomeView from './views/HomeView.vue'
import { defineComponent, h, ref } from 'vue'
import type { Component } from 'vue'
import { NSpace,NLayout,NLayoutSider,NMenu,NIcon,NLayoutFooter,NLayoutHeader,
  NScrollbar} from 'naive-ui'
import { RouterLink, RouterView } from 'vue-router'
import type { MenuOption } from 'naive-ui'
import {LaptopOutline,LogOutOutline,LogoElectron } from '@vicons/ionicons5'
import { Home,Tabler3DCubeSphere } from '@vicons/tabler'
import { Carbon,Workspace } from '@vicons/carbon'

function renderIcon (icon: Component) {
  return () => h(NIcon, null, { default: () => h(icon) })
}
const menuOptions: MenuOption[] = [
  {
    label: () =>h(RouterLink,{to: {name: 'home',}},{ default: () => 'Home' }),
    key: 'go-back-home',
    icon: renderIcon(Home)
  },
  {
    label: () =>h(RouterLink,{to: {name: 'ketcher',}},{ default: () => 'Ketcher' }),
    key: 'view-ketcher',
    icon: renderIcon(Carbon)
  },
  {
    label: () =>h(RouterLink,{to: {name: 'rdkit',}},{ default: () => 'RDKit' }),
    key: 'view-RDKit',
    icon: renderIcon(LogoElectron)
  },
  {
    label: () =>h(RouterLink,{to: {name: '3dmol',}},{ default: () => '3Dmol' }),
    key: 'view-3Dmol',
    icon: renderIcon(Tabler3DCubeSphere)
  },
  {
    label: () =>h(RouterLink,{to: {name: 'workspace',}},{ default: () => 'WorkSpace' }),
    key: 'view-workSpcae',
    icon: renderIcon(Workspace)
  }
]
</script>

<template>
<n-layout bordered position="absolute">
  <n-layout-header style="height: 60px; padding: 10px " bordered >
    <img class="logo" src="@/assets/logo.svg" width="40" height="40" />
  </n-layout-header>
  <n-layout has-sider position="absolute" style="top: 60px; bottom: 30px" bordered > 
        <n-layout-sider bordered
                        show-trigger
                        collapse-mode="width"
                        :collapsed-width="70"
                        :native-scrollbar="false"
                        :width="200">
          <n-menu :collapsed-width="70"
                  :collapsed-icon-size="28"
                  :options="menuOptions"/>
        </n-layout-sider>
        <n-layout content-style="padding: 20px" :native-scrollbar="false" bordered>
            <router-view>
            </router-view>
        </n-layout>
  </n-layout> 
  <n-layout-footer position="absolute" style="padding: 10px" bordered>
  </n-layout-footer>
</n-layout>
</template>

<style>
.logo {
  position:relative;
  left:10px;
}
</style>