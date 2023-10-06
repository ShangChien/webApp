<script setup lang="ts">
import { NStatistic } from 'naive-ui'
</script>

<template>
  <div class="flex flex-col items-center justify-center rd-2 bg-blue-50 h-90vh gap-4">
    <div class="w-40vw text-lg">
      <div class="i-fluent-emoji-thinking-face text-4xl inline-block" />
      <div class="inline">
        I'm rethinking the R&D style of molecule on the new web world
        just call me:
      </div>
      <div class="text-2xl inline">
        emm,
      </div>
      <div class="inline">
        which meaning: Exploring Molecular Magic, Efficient Molecular Mastery, Empowering Molecular Modeling...
      </div>
    </div>
    <div class="flex flex-col items-center justify-center bg-fuchsia-1 w-50vw rd-2">
      <div class="text-2xl ">
        Open Source Tools or Libraries used are as follows:
      </div>
      <div>
        <NStatistic label="平面分子画布">
          <a href="https://lifescience.opensource.epam.com/ketcher/index.html">Ketcher</a>
          (Apache License v2.0)
        </NStatistic>
      </div>
      <hr>
      <div>
        <NStatistic label="2D分子展示">
          <a href="https://greglandrum.github.io/rdkit-blog/">Rdkit</a>
          (BSD 3-Clause License)
        </NStatistic>
      </div>
      <hr>
      <div>
        <NStatistic label="3D分子展示">
          <a href="http://nglviewer.org/ngl/">NGL</a>
          (BSD-3-Clause license)
        </NStatistic>
      </div>
      <hr>
      <div>
        <NStatistic label="分子数据库">
          <a href="https://www.postgresql.org/about/">PostgreSQL</a>
          (Self-own License)
        </NStatistic>
      </div>
    </div>
  </div>
</template>

<style>
</style>
