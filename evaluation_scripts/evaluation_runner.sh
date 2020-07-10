#!/bin/bash

set -ex
graph_sets='small_graphs social_network web biological'
scenarios='full plateauBound withoutBucketQueue'

#for graph_set in $graph_sets; do
#  for scenario in $scenarios; do
#    python3 evaluation.py -g ${graph_set} -s ${scenario}|| true
#  done
#done

python3 evaluation.py -g small_graphs -s full
python3 evaluation.py -g small_graphs -s plateauBound
python3 evaluation.py -g small_graphs -s withoutBucketQueue
python3 evaluation.py -g social_network -s full
python3 evaluation.py -g social_network -s plateauBound
python3 evaluation.py -g social_network -s withoutBucketQueue
python3 evaluation.py -g web -s full
python3 evaluation.py -g web -s plateauBound
python3 evaluation.py -g biological -s full
python3 evaluation.py -g biological -s plateauBound
