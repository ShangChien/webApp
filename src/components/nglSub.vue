<script setup lang="ts">
import { ref,onMounted, onUpdated, onBeforeMount} from 'vue';
import  * as NGL from 'ngl/dist/ngl.js';
const viewport = ref(null);

function initViewer(viewport:any){
  let stage = new NGL.Stage(viewport.value,{backgroundColor:'white'});
  stage.loadFile("/moleculeData/1ycr.pdb")
    .then(function(o:any) {
      o.addRepresentation("cartoon");
      o.addRepresentation("licorice", { sele: "/0", multipleBond: "symmetric" })
      o.stage.setSpin(true);
      o.autoView();
      console.log(o);
      //悬停事件
      $(".nglView").hover(function(){
        o.stage.setSpin(false);
      },function(){
        o.stage.setSpin(true);
      });
     
      
      
  });
}
onMounted(()=>{
  initViewer(viewport)
})

</script>

<template>
<div ref="viewport" class="nglView"></div>
</template>

<style>
.nglView{
    width: 400px;
    height: 400px;
    position: relative;
}
</style>