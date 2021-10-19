#!/bin/bash

#This is for doing voxel calculation fo one tree: 1) the input tree (one tree)
#needs to be specified here, 2) this loops over sizes of voxels and number of parts
#Run identifier shown in the name of directory where results are stored is specified
#here  

numP=(1 3 5)
vox=(0.1 0.2 0.3 0.4) 
numPLen=${#numP[@]}
vLen=${#vox[@]}

folder=small
tree=Puut/this55-14.4267-14.2131-10.xml

echo Run folder ${folder} Tree ${tree}  NoPartsLevels ${numPLen}  Voxellevels ${vLen} 

  mkdir Resultapp/${folder}

    for (( l=0; l<${numPLen}; l++ ));
    do
      mkdir Resultapp/${folder}/numParts_${numP[$l]}
      echo Directory: Resultapp/${folder}/numParts_${numP[$l]}

      mkdir Resultapp/${folder}/numParts_${numP[$l]}/MeanStar
      echo Directory: Resultapp/${folder}/numParts_${numP[$l]}/MeanStar
      mkdir Resultapp/${folder}/numParts_${numP[$l]}/MeanStar/BoxDirNo
      echo Directory: Resultapp/${folder}/numParts_${numP[$l]}/MeanStar/BoxDirNo
      mkdir Resultapp/${folder}/numParts_${numP[$l]}/MeanStar/BoxDirYes
      echo Directory: Resultapp/${folder}/numParts_${numP[$l]}/MeanStar/BoxDirYes
      mkdir Resultapp/${folder}/numParts_${numP[$l]}/Directional
      echo Directory: Resultapp/${folder}/numParts_${numP[$l]}/Directional

      for (( k=0; k<${vLen};k++));
      do


     ./lig-radiation -self -inputTree ${tree} -resultfile Resultapp/${folder}/numParts_${numP[$l]}/MeanStar/BoxDirYes/boxdirYmean.dat -radMethod 3 -correctSTAR -boxDirEffect  -Voxboxside ${vox[$k]} -numParts ${numP[$l]} -appendMode


     ./lig-radiation -self -inputTree ${tree} -resultfile Resultapp/${folder}/numParts_${numP[$l]}/MeanStar/BoxDirNo/boxdirNmean.dat -radMethod 3  -correctSTAR  -Voxboxside ${vox[$k]} -numParts ${numP[$l]}  -appendMode


    ./lig-radiation -self -inputTree ${tree} -resultfile Resultapp/${folder}/numParts_${numP[$l]}/Directional/Boxdir.dat -radMethod 3  -correctSTAR  -calculateDirectionalStar -Voxboxside ${vox[$k]} -numParts ${numP[$l]}  -appendMode


   done
done

