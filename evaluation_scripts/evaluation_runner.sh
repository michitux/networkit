#!/bin/bash

set -ex
graph_sets='small_graphs social_network web biological'
scenarios='full plateauBound withoutBucketQueue'

#for graph_set in $graph_sets; do
#  for scenario in $scenarios; do
#    python3 evaluation.py -g ${graph_set} -s ${scenario}|| true
#  done
#done

python3 evaluation.py -g small_graphs -s full -o  || true
python3 evaluation.py -g small_graphs -s plateauBound -o || true
python3 evaluation.py -g small_graphs -s withoutBucketQueue -o || true

python3 evaluation.py -g biological -s full -o || true
python3 evaluation.py -g biological -s withoutBucketQueue -o || true

python3 evaluation.py -g web -s full -o || true

python3 evaluation.py -g social_network -s full -o || true

python3 evaluation.py -g biological -s plateauBound -o || true

