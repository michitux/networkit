import sys, getopt
import os
import networkit as nk
import timeit
import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser(prog='evaluation.py')
parser.add_argument('-g', '--graph_name')
parser.add_argument('-p', '--path')
parser.add_argument('-s', '--scenario', choices=['full', 'plateauBound', 'withoutBucketQueue', 'simple', 'init'])
parser.add_argument('-r', '--random_seed', type=int)
parser.add_argument('-o', '--overwrite', action='store_true')

args = vars(parser.parse_args())
graph_name = args['graph_name']
scenario = args['scenario']
overwrite = args['overwrite']
seed = args['random_seed']
input_path = "../input/"
output_path = args['path']
nk.setSeed(seed, False)
        
        
def getInitName(i):
    if (i == 0):
        return 'trivial'
    if (i == 1):
        return 'editing'
    if (i == 2):
        return 'random_insert'
    if (i == 3):
        return 'asc_degree_insert'
        
def executeMover (G, graph_name, init, m, s, r, p = 0, b_queue = True):
    mover = nk.community.QuasiThresholdEditingLocalMover(G, init, m, s, r, p, b_queue)
    a = timeit.default_timer()  
    mover.run()   
    delta = timeit.default_timer() - a
    edits = mover.getNumberOfEdits()
    usedIterations = mover.getUsedIterations()
    time = delta * 1000
    if(r):
        actualPlateau =  mover.getPlateauSize()
    else:
        actualPlateau = 0
    return [graph_name, G.numberOfNodes(), getInitName(init), m, s, r, p, 
            edits, usedIterations, actualPlateau, time]
            
def runOnGraph(graph_name, df):
    name = graph_name.split('/')[-1].split('.')[0]
    i = len(df.index)
    graph_path = input_path + graph_name
    if(graph_name.split('/')[0] ==  "facebook100"):
        G = nk.graphio.readMat(input_path + graph_name, key="A")
    if(graph_name.split('.')[1] == "graph"):
        G = nk.readGraph(graph_path, nk.Format.METIS)
    if(graph_name.split('.')[1] == "edgelist"):
        G = nk.readGraph(graph_path, nk.Format.SNAP)
    if(graph_name.split('.')[1] == "pairs"):
        G = nk.readGraph(graph_path, nk.Format.SNAP)
    G.indexEdges()    
    for init in initializations:
        for m in maxIterations:
            for s in sortPaths:
                for r in randomness:
                    if(r):
                        for p in plateauSize:
                            df.loc[i] = executeMover(G, name, init, m, s, r, p, b_queue)
                            i += 1
                    else:
                        df.loc[i] = executeMover(G, name, init, m, s, r, 0, b_queue)
                        i+=1
    return df



if(scenario == 'init'):
    maxIterations = [0]
    initializations = [0, 1, 2, 3]
    sortPaths = [True]
    randomness = [False]
    plateauSize = [0]
    b_queue = True
if(scenario == 'simple'):
    maxIterations = [5]
    initializations = [3]
    sortPaths = [True]
    randomness = [False]
    plateauSize = [0]
    b_queue = True
if(scenario == 'full'):
    maxIterations = [0, 5, 100]
    initializations = [0, 1, 2, 3]
    sortPaths = [True, False]
    randomness = [True, False]
    plateauSize = [5]
    b_queue = True
if(scenario == 'plateauBound'):
    initializations = [1, 3]
    maxIterations = [100]
    sortPaths = [True, False]
    randomness = [True]
    plateauSize = [1, 5, 100]
    b_queue = True
if(scenario == 'withoutBucketQueue'):
    initializations = [0, 1, 2, 3]
    maxIterations = [0, 5, 100]
    sortPaths = [True, False]
    randomness = [True, False]
    plateauSize = [5]
    b_queue = False
    
df = pd.DataFrame(columns  = ['graph', 
                              'n',
                              'initialization', 
                              'maxIterations', 
                              'sortPaths', 
                              'randomness', 
                              'plateauSize',
                              'edits',
                              'usedIterations',
                              'actualPlateau',
                              'time'])
                              
                              
graph_name_simple = graph_name.split('/')[-1].split('.')[0]
output_path += graph_name_simple + '/'

df = runOnGraph(graph_name, df)

if not os.path.exists(output_path):
    os.makedirs(output_path)
df['maxIterations'] = df['maxIterations'].apply(np.int64)
df['plateauSize'] = df['plateauSize'].apply(np.int64)
df['edits'] = df['edits'].apply(np.int64)
df['usedIterations'] = df['usedIterations'].apply(np.int64)
df['actualPlateau'] = df['actualPlateau'].apply(np.int64)
df['n'] = df['n'].apply(np.int64)
df.to_csv(output_path + graph_name_simple + '_' + scenario + '_' + str(seed) + '.csv', sep=',', encoding='utf-8')




