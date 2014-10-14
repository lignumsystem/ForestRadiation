#!/bin/bash
den=(23)
seed=(787237 536199  1676779  546327  235663 787231 536197 1576775 446327 135663)

dLen=${#den[@]}
sLen=${#seed[@]}

for (( i=0; i<${dLen}; i++ ));
  do
    mkdir Resultapp/${den[$i]}_accurate
    echo Directory: Resultapp/${den[$i]}_accurate
  for (( j=0; j<${sLen}; j++ ));
   do
      ./lig-radiation -generateLocations ${den[$i]} -manyTrees manytrees.txt -resultfile Resultapp/${den[$i]}_accurate/accurate_${seed[$j]}_den${den[$i]}.dat -seed ${seed[$j]} -radMethod 2  -treeDist 0.3 -X 20 -Y 20 -Z 10 -appendMode -GapRadius 0.6
  done
done

