#!/bin/bash

set -ex
graph_sets='small_graphs social_network web biological'
scenarios='full plateauBound withoutBucketQueue'
seeds='0 1 2 3 4 5 6 7 8 9'
limit=400
scenario=full
graph_set='biological'


small_graphs='adjnoun.graph dolphins.graph email.graph karate.graph lesmis.graph football.graph terrorist.edgelist grass_web.pairs'
#social_network='amazon.edgelist facebook100/Caltech36.mat dblp.edgelist youtube.edgelist facebook100/Penn94.mat lj.edgelist orkut.edgelist'
social_network1='facebook100/Caltech36.mat amazon.edgelist'
social_network2='dblp.edgelist facebook100/Penn94.mat youtube.edgelist'
#web='cnr-2000.graph in-2004.graph eu-2005.graph uk-2002.graph'
web='cnr-2000.graph in-2004.graph eu-2005.graph'

#python3 bionames.py -l ${limit}
#input="bio${limit}.txt"
#while IFS= read -r graph
#do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p "../output/QTM/biological${limit}/temp_${scenario}/" -s ${scenario} -r ${seed} &
#  done
#  wait
#done < "$input"
#python3 means.py -p "../output/QTM/biological${limit}/temp_${scenario}/"
#rm ${input}


graph_set='social_network'
scenario='simple'
for graph in $social_network1; do
  for seed in $seeds; do
    python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r ${seed} &
  done
  wait
done
python3 means.py -p "../output/QTM/${graph_set}/temp_${scenario}/"

    

graph_set='social_network'
scenario='full'
for graph in $social_network1; do
  for seed in $seeds; do
    python3 evaluation.py -g ${graph} -p "../output/QTM/${graph_set}/temp_${scenario}/" -s ${scenario} -r ${seed} &
  done
  wait
done
python3 means.py -p "../output/QTM/${graph_set}/temp_${scenario}/"




#for graph in $small_graphs; do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p '../output/QTM/small_graphs/temp_plateauBound/' -s 'plateauBound' -r ${seed} &
#  done
#  wait
#done
#python3 means.py -p "../output/QTM/small_graphs/temp_plateauBound/"

#for graph in $small_graphs; do
#  for seed in $seeds; do
#    python3 evaluation.py -g ${graph} -p '../output/QTM/small_graphs/temp_withoutBucketQueue/' -s 'withoutBucketQueue' -r ${seed} &
#  done
#  wait
#done
#python3 means.py -p "../output/QTM/small_graphs/temp_withoutBucketQueue/"


