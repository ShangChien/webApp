<script setup lang="ts">
import { NCol, NRow, NStatistic } from 'naive-ui'
</script>

<template>
  <div>
    <NRow>
      <NCol :span="20">
        <NStatistic label="平面分子画布">
          <a href="https://lifescience.opensource.epam.com/ketcher/index.html">Ketcher</a>
          (Apache License v2.0)
        </NStatistic>
      </NCol>
    </NRow>
    <hr>
    <NRow>
      <NCol :span="20">
        <NStatistic label="2D分子展示">
          <a href="https://greglandrum.github.io/rdkit-blog/">Rdkit</a>
          (BSD 3-Clause License)
        </NStatistic>
      </NCol>
    </NRow>
    <hr>
    <NRow>
      <NCol :span="20">
        <NStatistic label="3D分子展示">
          <a href="http://nglviewer.org/ngl/">NGL</a>
          (BSD-3-Clause license)
        </NStatistic>
      </NCol>
    </NRow>
    <hr>
    <NRow>
      <NCol :span="20">
        <NStatistic label="分子数据库">
          <a href="https://www.postgresql.org/about/">PostgreSQL</a>
          is released under the PostgreSQL License (similar to the BSD or MIT
          licenses)
        </NStatistic>
      </NCol>
    </NRow>
    <hr>
    <br>
    开源协议的主流类型:<br>
    <img src="/picture/openclass.png">
    <hr>
    <br>
    开源协议的通俗解释:<br>
    <img src="/picture/openbranch.png">
  </div>
</template>

<style>
</style>
