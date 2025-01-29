#!/bin/bash

for i in {1..40}; do
  
  filenameCPU="../runb382521/ana_${i}_tho_sao.nc"
  filenameGPU="ana_${i}_tho_sao.nc"
  
  if [[ -f $filenameCPU && -f $filenameGPU ]]; then
    echo "Compare files: $filenameCPU VS $filenameGPU:"
    cdo diffn "$filenameCPU" "$filenameGPU"
    echo # Empty string
  fi    

done

