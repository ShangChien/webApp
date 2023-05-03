import { createRouter, createWebHistory } from "vue-router";
import HomeView from "../views/HomeView.vue";

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: "/",
      name: "home",
      component: HomeView,
    },
    {
      path: "/ketcher",
      name: "ketcher",
      // route level code-splitting
      // this generates a separate chunk (About.[hash].js) for this route
      // which is lazy-loaded when the route is visited.
      component: () => import("../views/ketcher.vue"),
    },
    {
      path: "/molStore",
      name: "molStore",
      component: () => import("../views/molStore.vue"),
    },
    {
      path: "/3d",
      name: "3d",
      component: () => import("../views/threeDmol.vue"),
    },
    {
      path: "/task",
      name: "task",
      component: () => import("../views/task.vue"),
    },
    {
      path: "/enummole",
      name: "enumMole",
      component: () => import("../views/enumMole.vue"),
    },
    {
      path: "/enum",
      name: "enum",
      component: () => import("../views/enum.vue"),
    },
    {
      path: "/python",
      name: "python",
      component: () => import("../views/pythonView.vue"),
    },
  ],
});

export default router;
