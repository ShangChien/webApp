<script setup lang="ts">
import { ref, onMounted, watch, watchEffect } from 'vue';
import * as $3Dmol from '3dmol/build/3Dmol-nojquery-min.js';
const threeDmol = ref(null);

onMounted(()=>{
  $(function() {
    let element = threeDmol.value
    let config = { backgroundColor: 'white' };
    let viewer = $3Dmol.createViewer( element, config );
    let pdbUri = '/moleculeData/1ycr.pdb';
    jQuery.ajax( pdbUri, { 
      success: function(data) {
        let v = viewer;
        v.addModel( data, "pdb" );                       /* load data */
        v.setStyle({}, {cartoon: {color: 'spectrum'}});  /* style all atoms */
        v.zoomTo();                                      /* set camera */
        v.render();                                      /* render scene */
        v.spin('y',0.2);
        //悬停事件
        $(".mol-container").hover(function(){
          v.spin('y',0);
        },function(){
          v.spin('y',0.2);
        });
      },
      error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
      },
    });
  });
})


</script>
<template>
<div ref="threeDmol" class="mol-container"></div>
</template>
<style>
.mol-container {
  width: 100%;
  height: 400px;
  position: relative;
}
</style>