#!/bin/bash
vector_size=$1

rm -f hmmprototype.txt
mean=""
variance=""

         i="0"
         while [  $i -lt $vector_size ]; do
             mean="${mean}0.0 "
             variance="${variance}1.0 "
             i=$((i+1))
         done


echo "<BeginHMM>" >> hmmprototype.txt
cat <<EOT>> hmmprototype.txt
 <NumStates> 6 <VecSize> $vector_size <USER> <nullD> <diagC>
 <State> 2 <NumMixes> 1   
  <Mixture> 1 1.0
   <Mean> $vector_size
     $mean
   <Variance> $vector_size
     $variance
 <State> 3 <NumMixes> 1   
  <Mixture> 1 1.0
    <Mean> $vector_size
      $mean
    <Variance> $vector_size
      $variance
 <State> 4 <NumMixes> 1   
  <Mixture> 1 1.0
    <Mean> $vector_size
      $mean
    <Variance> $vector_size
      $variance
 <State> 5 <NumMixes> 1   
  <Mixture> 1 1.0
    <Mean> $vector_size
      $mean
    <Variance> $vector_size
      $variance
 <TransP> 6
  0.000e+0    1.000e+0    0.000e+0    0.000e+0    0.000e+0   0.000e+0
  0.000e+0    6.000e-1    4.000e-1    0.000e+0    0.000e+0   0.000e+0
  0.000e+0    0.000e+0    6.000e-1    4.000e-1    0.000e+0   0.000e+0
  0.000e+0    0.000e+0    0.000e+0    6.000e-1    4.000e-1   0.000e+0 
  0.000e+0    0.000e+0    0.000e+0    0.000e+0    9.000e-1   1.000e-1  
  0.000e+0    0.000e+0    0.000e+0    0.000e+0    0.000e+0   0.000e+1
<EndHMM>
EOT
