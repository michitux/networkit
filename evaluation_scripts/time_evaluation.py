import sys, getopt
import os
import networkit as nk
import timeit
import pandas as pd
import numpy as np
import argparse

#python3 time_evaluation.py -i 3 -m 100 -s True -r True -p 5 -b True

parser = argparse.ArgumentParser(prog='evaluation.py')
parser.add_argument('-i', '--init', type=int, choices=[0, 1, 2, 3])
parser.add_argument('-m', '--maxIterations', type=int)
parser.add_argument('-s', '--sortPaths', type=bool)
parser.add_argument('-r', '--randomness', type=bool)
parser.add_argument('-p', '--plateau', type=int)
parser.add_argument('-b', '--bucketQueue', type=bool)

args = vars(parser.parse_args())
init = args['init']
maxIterations = args['maxIterations']
sortPaths = args['sortPaths']
randomness = args['randomness']
plateau = args['plateau']
bucketQueue = args['bucketQueue']


graph_names = ['adjnoun.graph', 'dolphins.graph', 'karate.graph', 'lesmis.graph', 'football.graph', 'terrorist.edgelist', 'grass_web.pairs',
'jazz.graph', 'polbooks.graph', 'amazon.edgelist', 'facebook100/Caltech36.mat', 'dblp.edgelist', 'youtube.edgelist', 'facebook100/Penn94.mat',
'lj.edgelist', 'orkut.edgelist', 'cnr-2000.graph', 'in-2004.graph', 'eu-2005.graph', 'uk-2002.graph']
#graph_names = ['adjnoun.graph', 'dolphins.graph', 'karate.graph', 'lesmis.graph']
input_path = '../input/'
df = pd.DataFrame(columns  = ['graph', 'iteration', 'edits', 'time'])
row = 0
for graph_name in graph_names:
    name = graph_name.split('/')[-1].split('.')[0]
    i = len(df.index)
    graph_path = input_path + graph_name
    if(graph_name.split('/')[0] ==  "facebook100"):
        G = nk.graphio.readMat(input_path + graph_name, key="A")
    if(graph_name.split('.')[-1] == "graph"):
        G = nk.readGraph(graph_path, nk.Format.METIS)
    if(graph_name.split('.')[-1] == "edgelist"):
        G = nk.readGraph(graph_path, nk.Format.SNAP, continuous=False, directed=False)
    if(graph_name.split('.')[-1] == "pairs"):
        G = nk.readGraph(graph_path, nk.Format.SNAP)
    G.indexEdges()
    mover = nk.community.QuasiThresholdEditingLocalMover(G, init, maxIterations, sortPaths, randomness, plateau, bucketQueue)
    mover.run()
    edits = mover.getRunningInfo()[b'edits']
    time = mover.getRunningInfo()[b'time']
    for i in range(len(edits)):
        df.loc[row] = [graph_name, i, edits[i], time[i]]
        row += 1
print(df)
df['iteration'] = df['iteration'].apply(np.int64)
df['edits'] = df['edits'].apply(np.int64)
df['time'] = df['time'].apply(np.int64)

output_path = '../output/QTM/'
if not os.path.exists(output_path):
    os.makedirs(output_path)
file_name = output_path + 'time_' + str(init) + '_' + str(maxIterations) + '_' + str(sortPaths) + '_' + str(randomness) + '_' + str(plateau) + '_' + str(bucketQueue) + '.csv'
df.to_csv(file_name, sep=',', encoding='utf-8')
