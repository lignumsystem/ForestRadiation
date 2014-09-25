#!/bin/bash
# To run this script there is a need to have a folder named Resultapp that will then store all the other results in the Resultapp folder. If your choice of the name is differnt then please change it accordinly in all the places here. i.e. Find and Replace and it will do the job

den=(150) # This is where the densities of the tress are given
seed=(787237  536199  1676779  546327  235663 787231 536197 1576775 446327 135663) 
numP=(3)
vox=(0.1 0.2 0.3 0.4)
dLen=${#den[@]} #Length of the density
sLen=${#seed[@]} # Length of the seeds list
numPLen=${#numP[@]}
vLen=${#vox[@]} 
for (( i=0; i<${dLen}; i++ ));
 do
   mkdir Resultapp/${den[$i]}
    for (( l=0; l<${numPLen}; l++ ));
    do
      mkdir Resultapp/${den[$i]}/numParts_${numP[$l]}
   for (( j=0; j<${sLen}; j++ ));
   do
   echo ${den[$i]}
   echo ${seed[$j]}
     for (( k=0; k<${vLen};k++));
       do
     
    mkdir Resultapp/${den[$i]}/numParts_${numP[$l]}/MeanStar # make changes to the path where you need to store the data here
    mkdir Resultapp/${den[$i]}/numParts_${numP[$l]}/MeanStar/BoxDirYes  # Make changes to the path if needed here.
     echo ${vox[$k]}
     echo -----------------------------------------
     ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/numParts_${numP[$l]}/MeanStar/BoxDirYes/boxdirYmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR -boxDirEffect  -Voxboxside ${vox[$k]} -numParts ${numP[$l]} -zeroWoodyRadius -appendMode
 cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDY_${vox[$k]}.dat

    

#*********************************************************************************************************************************
    mkdir Resultapp/${den[$i]}/numParts_${numP[$l]}/MeanStar/BoxDirNo
   ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/numParts_${numP[$l]}/MeanStar/BoxDirNo/boxdirNmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR  -Voxboxside ${vox[$k]} -numParts ${numP[$l]}  -appendMode
  cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDN_${vox[$k]}.dat

    
#*********************************************************************************************************************************
   mkdir Resultapp/${den[$i]}/numParts_${numP[$l]}/Directional # Make changes to the path if needed here.
    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/numParts_${numP[$l]}/Directional/boxdirNdir_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR  -calculateDirectionalStar -Voxboxside ${vox[$k]} -numParts ${numP[$l]}  -appendMode
 cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_DIR_${vox[$k]}.dat

      done
    done
   done
done

