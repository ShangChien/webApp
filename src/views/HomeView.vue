<script setup lang="ts">
import { NCol,NRow,NStatistic } from 'naive-ui'
</script>

<template>
  <n-row>
    <n-col :span="20">
      <n-statistic label="画布" >
        <a href='https://lifescience.opensource.epam.com/ketcher/index.html'>Ketcher</a>
        (Apache License v2.0)
      </n-statistic>
    </n-col>
  </n-row>
  <hr>
  <n-row>
    <n-col :span="20">
      <n-statistic label="2D.js" >
        <a href='https://greglandrum.github.io/rdkit-blog/'>Rdkit</a>
        (BSD 3-Clause License)
      </n-statistic>
    </n-col>
  </n-row>
  <hr>
  <n-row>
    <n-col :span="20">
      <n-statistic label="3D.js" >
        <a href='http://3dmol.csb.pitt.edu/index.html'>3Dmol</a>
        (BSD-3-Clause license)
      </n-statistic>
    </n-col>
  </n-row>
  <hr>
  <n-row>
    <n-col :span="20">
      <n-statistic label="Data-base" >
        <a href='https://www.postgresql.org/about/'>PostgreSQL</a>
        is released under the PostgreSQL License (similar to the BSD or MIT licenses)
      </n-statistic>
    </n-col>
  </n-row>
  <hr>
  <br>
  开源协议的主流类型:<br>
  <img src="/picture/openclass.png"  />
  <hr>
  <br>
  开源协议的通俗解释:<br>
  <img src="/picture/openbranch.png"  />
</template>

