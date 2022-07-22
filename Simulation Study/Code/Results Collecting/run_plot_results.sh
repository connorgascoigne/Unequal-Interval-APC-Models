#!/bin/bash    

for family in {1..3} # gaussian, poisson, binomial
do

  for dataMod in {1..4} # apc, ap, ac, a
  do
  
    for aggregation in {1..3} # M = 1, 3 and 5
    do
    
      echo "$family" "$dataMod" "$aggregation"
      Rscript --vanilla plot_results.R $family $dataMod $aggregation
    
    done
  
  done

done