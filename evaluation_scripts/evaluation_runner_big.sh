#!/bin/bash

set -ex
graph_sets='social_network web'
scenarios='full plateauBound'

declare -A graphs
graphs=(["social_network"]='lj.edgelist orkut.edgelist'
["web"]='uk-2002.graph')


for graph_set in $graph_sets; do
  for scenario in $scenarios; do
    graph_list=${graphs[$graph_set]}
    for graph in $graph_list; do
      taskset -c 0 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 0 &
      taskset -c 1 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 1 &
      taskset -c 2 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 2 &
      taskset -c 3 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 3 &
      taskset -c 4 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 4 &
      taskset -c 8 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 5 &
      taskset -c 9 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 6 &
      taskset -c 10 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 7 &
      taskset -c 11 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 8 &
      taskset -c 12 python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r 9 &
      wait
    done
    python3 means.py -p "../output/QTM/${graph_set}/temp_${scenario}/"
  done
done



    




