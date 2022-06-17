<script setup lang="ts">
import {
  NSpace,
  NTable,
  NThead,
  NTr,
  NTh,
  NTbody,
  NTd,
  NProgress,
  NButtonGroup,
  NButton,
  NIcon,
} from "naive-ui";
import { Delete, View, Edit } from "@vicons/carbon";
</script>

<template>
  <n-space vertical>
    <n-table :single-line="false" size="small">
      <n-thead>
        <n-tr>
          <n-th>序号</n-th>
          <n-th>名称</n-th>
          <n-th>类型</n-th>
          <n-th>状态</n-th>
          <n-th>操作</n-th>
        </n-tr>
      </n-thead>
      <n-tbody>
        <n-tr>
          <n-td>1</n-td>
          <n-td>220510-ETA-0028</n-td>
          <n-td>结构优化</n-td>
          <n-td>
            <n-progress
              type="line"
              :percentage="60"
              :indicator-placement="'inside'"
              color="#52c234"
              processing
            />
          </n-td>
          <n-td>
            <n-button-group>
              <n-button size="small" round tertiary type="primary">
                <template #icon>
                  <n-icon><View /></n-icon>
                </template>
                查看
              </n-button>
              <n-button size="small" round tertiary type="info">
                <template #icon>
                  <n-icon><Edit /></n-icon>
                </template>
                编辑
              </n-button>
              <n-button size="small" round tertiary type="error">
                <template #icon>
                  <n-icon><Delete /></n-icon>
                </template>
                删除
              </n-button>
            </n-button-group>
          </n-td>
        </n-tr>
        <n-tr>
          <n-td>...</n-td>
          <n-td>...</n-td>
          <n-td>...</n-td>
          <n-td>...</n-td>
          <n-td>...</n-td>
        </n-tr>
        <n-tr>
          <n-td>...</n-td>
          <n-td>...</n-td>
          <n-td>...</n-td>
          <n-td>...</n-td>
          <n-td>...</n-td>
        </n-tr>
      </n-tbody>
    </n-table>
  </n-space>
</template>
