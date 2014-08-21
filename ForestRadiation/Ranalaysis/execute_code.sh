#!/bin/bash
# To run this script there is a need to have a folder named Resultapp that will then store all the other results in the Resultapp folder. If your choice of the name is differnt then please change it accordinly in all the places here. i.e. Find and Replace and it will do the job

den=(150) # This is where the densities of the tress are given
seed=(787237 536199  1676779  546327  235663 787231 536197 1576775 446327 135663) # seeds are specified
Voxsize=(0.1 0.2 0.3 0.4)
dLen=${#den[@]} #Length of the density
sLen=${#seed[@]} # Length of the seeds list
vLen=${#Voxsize[@]} # Length of the Voxsize list
for (( i=0; i<${dLen}; i++ ));
 do
   mkdir Resultapp/${den[$i]}
   for (( j=0; j<${sLen}; j++ ));
   do
     for((k = 0;k<${vLen};k++));
   do 
   echo ${den[$i]}
   echo ${seed[$j]}
    
    mkdir Resultapp/${den[$i]}/MeanStar # make changes to the path where you need to store the data here
    mkdir Resultapp/${den[$i]}/MeanStar/BoxDirYes  # Make changes to the path if needed here.
     ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirYes/boxdirYmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR -boxDirEffect  -Voxboxside ${Voxsize[$k]}  -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDY_${Voxsize[$k]} .dat
#*********************************************************************************************************************************
    mkdir Resultapp/${den[$i]}/MeanStar/BoxDirNo
        ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirNo/boxdirNmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR  -Voxboxside${Voxsize[$k]}  -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDN_${Voxsize[$k]}.dat
#*********************************************************************************************************************************
   mkdir Resultapp/${den[$i]}/Directional # Make changes to the path if needed here.
     ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/Directional/boxdirNdir_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR  -calculateDirectionalStar -Voxboxside ${Voxsize[$k]} -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_DIR_${Voxsize[$k]}.dat

    done
   done
done
