<script setup lang="ts">
import { h } from 'vue'
import type { Component } from 'vue'
import { NLayout,NLayoutSider,NMenu,NIcon,NLayoutFooter,NLayoutHeader,NMessageProvider} from 'naive-ui'
import { RouterLink, RouterView } from 'vue-router'
import type { MenuOption } from 'naive-ui'
import { LogoElectron } from '@vicons/ionicons5'
import { Home,Tabler3DCubeSphere,BrandSlack } from '@vicons/tabler'
import { Carbon,Data1 } from '@vicons/carbon'


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
    label: () =>h(RouterLink,{to: {name: '3d',}},{ default: () => '3D' }),
    key: 'view-3D',
    icon: renderIcon(Tabler3DCubeSphere)
  },
  {
    label: () =>h(RouterLink,{to: {name: 'task',}},{ default: () => 'Task' }),
    key: 'view-task',
    icon: renderIcon(Data1)
  },
  {
    label: () =>h(RouterLink,{to: {name: 'enumMole',}},{ default: () => 'Enum Molecule' }),
    key: 'view-workSpcae',
    icon: renderIcon(BrandSlack)
  }
]
</script>

<template>
<n-layout bordered position="absolute" :native-scrollbar="false" >
  <n-layout-header style="height: 60px; padding: 1px " bordered >
    <img class="logo" src="/picture/colorlogo.svg" width="55"  />
  </n-layout-header>
  <n-layout has-sider position="absolute" style="top: 60px" bordered > 
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
            <n-message-provider>
              <router-view v-slot="{ Component }">
                <keep-alive>
                  <component :is="Component" />
                </keep-alive>
              </router-view>
            </n-message-provider>
        </n-layout>
  </n-layout> 
</n-layout>
</template>

<style>
.logo{
  position: relative;
  left: 10px;
}
</style>