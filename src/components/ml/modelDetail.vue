<script setup lang='ts'>
import { NCheckbox, NCheckboxGroup, NCollapse, NCollapseItem, NRadio, NRadioGroup, NSpace } from 'naive-ui'

const task = defineModel<string>('task')
const taskList = ['regression', 'classification', 'multiclass', 'multilabel_classification', 'multilabel_regression']

const models = defineModel<string[]>('models')
const modelsList = ['HOMO', 'LUMO', 'Eg', 'Representation']

const fileType = defineModel<string>('fileType')
const fileTypeList = ['*.mol(3D)', '*.mol(2D)', '*.smi']
</script>

<template>
  <NCollapse class="noCopy" default-expanded-names="1" accordion>
    <template #arrow>
      <div class="i-icon-park-config text-3xl box-border" />
    </template>
    <NCollapseItem class="box-border" name="1">
      <template #header>
        <div class="text-xl box-border">
          Setting
        </div>
      </template>
      <div class="flex flex-col justify-around items-center gap-1 box-border">
        <div class="flex-auto flex justify-start items-center gap-3 w-full box-border">
          <div class="w-100px text-right flex-none">
            task:
          </div>
          <NRadioGroup v-model:value="task" name="radiogroup" class="bg-sky-1 px-2 py-1 rd-2 box-border">
            <NSpace>
              <NRadio
                v-for="(i_task, index) in taskList"
                :key="index" :value="i_task" class="p-1"
                :disabled="i_task !== 'regression'"
                :class="[i_task === task ? 'bg-blue-3 rd-1 ' : '']"
              >
                {{ i_task }}
              </NRadio>
            </NSpace>
          </NRadioGroup>
        </div>
        <div class="flex-auto flex justify-start items-center gap-3 w-full box-border">
          <div class="w-100px text-right flex-none">
            models:
          </div>
          <NCheckboxGroup v-model:value="models" class="bg-sky-1 p-1 rd-2 box-border">
            <NCheckbox
              v-for="(i_model, index) in modelsList"
              :key="index" :value="i_model"
              :disabled="i_model === 'Representation'"
              :label="i_model"
              class="mx-1 pl-2 p-1"
              :class="[models.includes(i_model) ? 'bg-blue-3 rd-1 ' : '']"
            />
          </NCheckboxGroup>
          <div class="text-nowrap">
            (multi select)
          </div>
        </div>
        <div class="flex-auto flex justify-start items-center gap-3 w-full box-border">
          <div class="w-100px text-right flex-none">
            fileType:
          </div>
          <NRadioGroup v-model:value="fileType" name="fileType" class="bg-sky-1 px-2 py-1 rd-2 box-border">
            <NSpace>
              <NRadio v-for="(i_fileType, index) in fileTypeList" :key="index" :value="i_fileType" class="p-1" :class="[i_fileType === fileType ? 'bg-blue-3 rd-1 ' : '']">
                {{ i_fileType }}
              </NRadio>
            </NSpace>
          </NRadioGroup>
        </div>
      </div>
    </NCollapseItem>
  </NCollapse>
</template>

<style scoped>
.noCopy {
-webkit-touch-callout:none;  /*系统默认菜单被禁用*/
-webkit-user-select:none; /*webkit浏览器*/
-khtml-user-select:none; /*早期浏览器*/
-moz-user-select:none;/*火狐*/
-ms-user-select:none; /*IE10*/
user-select:none;
}
</style>
