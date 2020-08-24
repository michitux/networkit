#!/bin/bash

set -ex
graph_sets='small_graphs'
#scenarios='full plateauBound withoutBucketQueue'
scenarios='full plateauBound'
seeds='0 1 2 3 4 5 6 7 8 9'

declare -A graphs
graphs=(["small_graphs"]='adjnoun.graph dolphins.graph karate.graph lesmis.graph football.graph terrorist.edgelist grassweb.pairs jazz.graph polbooks.graph'
["social_network"]='amazon.edgelist facebook100/Caltech36.mat dblp.edgelist youtube.edgelist facebook100/Penn94.mat'
["web"]='cnr-2000.graph in-2004.graph eu-2005.graph')

#small_graphs='adjnoun.graph dolphins.graph karate.graph lesmis.graph football.graph terrorist.edgelist grass_web.pairs jazz.graph polbooks.graph'
#social_network='amazon.edgelist facebook100/Caltech36.mat dblp.edgelist youtube.edgelist facebook100/Penn94.mat lj.edgelist orkut.edgelist'
#web='cnr-2000.graph in-2004.graph eu-2005.graph uk-2002.graph'

#graph_set='generated'
#input="generatednames.txt"
#scenario='full'
#while IFS= read -r graph
#do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p "../output/QTM_3/${graph_set}/temp_${scenario}/" -s ${scenario} -r ${seed} &
#  done
#  wait
#done < "$input"
#python3 means.py -p "../output/QTM_3/${graph_set}/temp_${scenario}/"

#scenario='plateauBound'
#while IFS= read -r graph
#do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p "../output/QTM_3/${graph_set}/temp_${scenario}/" -s ${scenario} -r ${seed} &
#  done
#  wait
#done < "$input"
#python3 means.py -p "../output/QTM_3/${graph_set}/temp_${scenario}/"


for graph_set in $graph_sets; do
  for scenario in $scenarios; do
    graph_list=${graphs[$graph_set]}
    for graph in $graph_list; do
      for seed in $seeds; do
        python3 evaluation.py -g ${graph} -p "../output/QTM_3/${graph_set}/temp_${scenario}/" -s ${scenario} -r ${seed} &
      done
      wait
    done
    python3 means.py -p "../output/QTM_3/${graph_set}/temp_${scenario}/"
  done
done

#limit=400
#input="bio${limit}.txt"
#scenario='full'
#while IFS= read -r graph
#do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p "../output/QTM_3/biological${limit}/temp_${scenario}/" -s ${scenario} -r ${seed} &
#  done
#  wait
#done < "$input"
#python3 means.py -p "../output/QTM_3/biological${limit}/temp_${scenario}/"

#scenario='plateauBound'
#while IFS= read -r graph
#do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p "../output/QTM_3/biological${limit}/temp_${scenario}/" -s ${scenario} -r ${seed} &
#  done
#  wait
#done < "$input"
#python3 means.py -p "../output/QTM_3/biological${limit}/temp_${scenario}/"
