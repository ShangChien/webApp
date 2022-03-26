import { createRouter, createWebHistory } from 'vue-router'
import HomeView from '../views/HomeView.vue'

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: '/',
      name: 'home',
      component: HomeView
    },
    {
      path: '/ketcher',
      name: 'ketcher',
      // route level code-splitting
      // this generates a separate chunk (About.[hash].js) for this route
      // which is lazy-loaded when the route is visited.
      component: ()=>(import('../views/ketcher.vue'))
    },
    {
      path: '/rdkit',
      name: 'rdkit',
      component: ()=>(import('../views/rdkit.vue'))
    },
    {
      path: '/3dmol',
      name: '3dmol',
      component: ()=>(import('../views/threeDmol.vue')),
    },
    {
      path: '/workspace',
      name: 'workspace',
      component: ()=>(import('../views/workSpace.vue'))
    },
  ]
})

export default router
