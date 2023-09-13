<script setup lang="ts">
import { h, ref } from 'vue'
import type { Component } from 'vue'
import {
  NIcon,
  NLayout,
  NLayoutFooter,
  NLayoutHeader,
  NLayoutSider,
  NMenu,
  NMessageProvider,
} from 'naive-ui'
import { RouterLink, RouterView } from 'vue-router'
import { LogoElectron } from '@vicons/ionicons5'

// eslint-disable-next-line unused-imports/no-unused-imports
import { BrandSlack, Home, Tabler3DCubeSphere } from '@vicons/tabler'
import { Carbon, Data1, DataVis1 } from '@vicons/carbon'
import ReloadPrompt from '@/components/ReloadPrompt.vue'

const collapsed = ref(true)
function renderIcon(icon: Component) {
  return () => h(NIcon, null, { default: () => h(icon) })
}
const menuOptions = [
  {
    label: () =>
      h(RouterLink, { to: { name: 'home' } }, { default: () => 'Home' }),
    key: 'go-back-home',
    icon: renderIcon(Home),
  },
  {
    label: () =>
      h(RouterLink, { to: { name: 'ketcher' } }, { default: () => 'Ketcher' }),
    key: 'view-ketcher',
    icon: renderIcon(Carbon),
  },
  {
    label: () =>
      h(RouterLink, { to: { name: 'molStore' } }, { default: () => 'molStore' }),
    key: 'view-molStore',
    icon: renderIcon(LogoElectron),
  },
  {
    label: () => h(RouterLink, { to: { name: '3d' } }, { default: () => '3D' }),
    key: 'view-3D',
    icon: renderIcon(Tabler3DCubeSphere),
  },
  {
    label: () =>
      h(RouterLink, { to: { name: 'task' } }, { default: () => 'Task' }),
    key: 'view-task',
    icon: renderIcon(Data1),
  },
  // {
  //   label: () =>
  //     h(RouterLink, { to: { name: 'enumMole' } }, { default: () => 'Enum Molecule' }),
  //   key: 'view-workSpcae',
  //   icon: renderIcon(BrandSlack),
  // },
  {
    label: () =>
      h(RouterLink, { to: { name: 'enum' } }, { default: () => 'Enum' }),
    key: 'view-enum',
    icon: renderIcon(DataVis1),
  },
  {
    label: () =>
      h(RouterLink, { to: { name: 'python' } }, { default: () => 'python' }),
    key: 'view-python',
    icon: () => h(NIcon, null, { default: () => h('div', { class: 'i-tabler-brand-python' }) }),
  },
  {
    label: () =>
      h(RouterLink, { to: { name: 'ml' } }, { default: () => 'ml' }),
    key: 'view-ml',
    icon: () => h('img', { class: 'h-full', alt: 'mechine learning', src: 'src/assets/core-ml-256.png' }),
  },
]
</script>

<template>
  <NLayout class="font-lx-b" bordered position="absolute" style="height: 100vh" :native-scrollbar="false" @contextmenu.prevent>
    <NLayoutHeader style="height: 6vh; padding: 2px" bordered>
      <img class="logo w-5vh" src="/favicon.ico">
      <ReloadPrompt />
    </NLayoutHeader>
    <NLayout has-sider position="absolute" style="top: 6vh; bottom: 2vh" bordered>
      <NLayoutSider
        bordered
        collapse-mode="width"
        :collapsed="collapsed"
        :collapsed-width="70"
        :native-scrollbar="false"
        :width="200"
      >
        <NMenu
          :collapsed-width="70"
          :collapsed-icon-size="28"
          :options="menuOptions"
          :collapsed="collapsed"
        />
        <div class="w-100% text-center position-absolute bottom-0">
          <div
            class="rd-1 ma-2 hover:(bg-violet-100 text-violet cursor-pointer duration-210 ease-in-out)"
            @click="collapsed = !collapsed"
          >
            <div v-if="collapsed" class="i-tabler-layout-sidebar-left-expand text-3xl ma-1" />
            <div v-else class="i-tabler-layout-sidebar-right-expand text-3xl ma-1" />
          </div>
        </div>
      </NLayoutSider>
      <NLayout
        class="p-1vh pr-0 box-border"
        :native-scrollbar="false"
        bordered
      >
        <NMessageProvider>
          <RouterView v-slot="{ Component: comp }" class="mr-3 ">
            <Transition name="fade" mode="out-in">
              <keep-alive>
                <component :is="comp" />
              </keep-alive>
            </Transition>
          </RouterView>
        </NMessageProvider>
      </NLayout>
    </NLayout>
    <NLayoutFooter
      bordered
      position="absolute"
      style="height: 2vh; padding: 2px"
    />
  </NLayout>
</template>

<style>
.logo {
  position: relative;
  left: 10px;
}

.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.2s ease;
}

.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}
.logo {
-webkit-touch-callout:none;  /*系统默认菜单被禁用*/
-webkit-user-select:none; /*webkit浏览器*/
-khtml-user-select:none; /*早期浏览器*/
-moz-user-select:none;/*火狐*/
-ms-user-select:none; /*IE10*/
user-select:none;
}
</style>
