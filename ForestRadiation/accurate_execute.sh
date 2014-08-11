#!/bin/bash
den=(300)
seed=(787237  536199  1676779  546327  235663)
dLen=${#den[@]}
sLen=${#seed[@]}

for (( i=0; i<${dLen}; i++ ));
  do
    mkdir Resultapp/${den[$i]} 
  for (( j=0; j<${sLen}; j++ )); 
   do 
      ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}/accurate_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 2 -zeroWoodyRadius -treeDist 0.3 -X 20 -Y 20 -Z 10 -appendMode
   done
done

