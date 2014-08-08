#!/bin/bash
# To run this script there is a need to have a folder named Resultapp that will then store all the other results in the Resultapp folder. If your choice of the name is differnt then please change it accordinly in all the places here. i.e. Find and Replace and it will do the job

den=( 100 300) # This is where the densities of the tress are given
seed=(787237 536199  1676779  546327  235663) # seeds are specified
dLen=${#den[@]} #Length of the density
sLen=${#seed[@]} # Length of the seeds list
for (( i=0; i<${dLen}; i++ ));
 do
   mkdir Resultapp/${den[$i]}
   for (( j=0; j<${sLen}; j++ ));
   do
   echo ${den[$i]}
   echo ${seed[$j]}
    
    mkdir Resultapp/${den[$i]}/MeanStar # make changes to the path where you need to store the data here
    mkdir Resultapp/${den[$i]}/MeanStar/BoxDirYes  # Make changes to the path if needed here.
     ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirYes/boxdirYmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR -boxDirEffect  -Voxboxside 0.1 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDY_0.1.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirYes/boxdirYmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10   -correctSTAR -boxDirEffect  -Voxboxside 0.2 -zeroWoodyRadius -appendMode
cp runfile.dat runfile_${seed[$j]}_${den[$i]}_BDY_0.2.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirYes/boxdirYmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10  -correctSTAR -boxDirEffect  -Voxboxside 0.3 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDY_0.3.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirYes/boxdirYmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3 -treeDist 0.3 -X 20 -Y 20 -Z 10   -correctSTAR -boxDirEffect  -Voxboxside 0.4 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDY_0.4.dat

#*********************************************************************************************************************************
    mkdir Resultapp/${den[$i]}/MeanStar/BoxDirNo
   ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirNo/boxdirNmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR  -Voxboxside 0.1 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDN_0.1.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirNo/boxdirNmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10   -correctSTAR  -Voxboxside 0.2 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDN_0.2.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirNo/boxdirNmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10  -correctSTAR  -Voxboxside 0.3 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDN_0.3.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/MeanStar/BoxDirNo/boxdirNmean_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10   -correctSTAR  -Voxboxside 0.4 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_BDN_0.4.dat

#*********************************************************************************************************************************
   mkdir Resultapp/${den[$i]}/Directional # Make changes to the path if needed here.
    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/Directional/boxdirNdir_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10 -correctSTAR  -calculateDirectionalStar -Voxboxside 0.1 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_DIR_0.1.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/Directional/boxdirNdir_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10   -correctSTAR   -calculateDirectionalStar -Voxboxside 0.2 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_DIR_0.2.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/Directional/boxdirNdir_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3  -treeDist 0.3 -X 20 -Y 20 -Z 10  -correctSTAR -calculateDirectionalStar -Voxboxside 0.3 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_DIR_0.3.dat

    ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/Directional/boxdirNdir_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 3 -treeDist 0.3 -X 20 -Y 20 -Z 10   -correctSTAR  -calculateDirectionalStar -Voxboxside 0.4 -zeroWoodyRadius -appendMode
cp runfile.dat Resultapp/${den[$i]}/runfile_${seed[$j]}_${den[$i]}_DIR_0.4.dat
   done
done

