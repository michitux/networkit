import sys, getopt
import os
import networkit as nk
import timeit
import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser(prog='evaluation.py')
parser.add_argument('-g', '--graph_set', choices=['small_graphs', 'biological', 'social_network', 'web', 'generated'])
parser.add_argument('-s', '--scenario', choices=['full', 'plateauBound', 'withoutBucketQueue'])
parser.add_argument('-o', '--overwrite', action='store_true')

args = vars(parser.parse_args())
print(args)
graph_set = args['graph_set']
scenario = args['scenario']
overwrite = args['overwrite']
seeds = list(range(9))
input_path = "../input/"
output_path = '../output/QTG/' + graph_set + '/' + scenario + '/'
        
        
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
    meanActualPlateau = 0.0
    meanEdits = 0.0
    meanTime = 0.0
    meanUsedIterations = 0.0
    minEdits = G.upperNodeIdBound()*G.upperNodeIdBound()
    maxUsedIterations = 0
    maxActualPlateau = 0
    for seed in seeds:
        nk.setSeed(seed, False)
        mover = nk.community.QuasiThresholdEditingLocalMover(G, init, m, s, r, p, b_queue)
        a = timeit.default_timer()  
        mover.run()   
        delta = timeit.default_timer() - a
        meanEdits = meanEdits + mover.getNumberOfEdits()
        minEdits = min(minEdits, mover.getNumberOfEdits())
        meanUsedIterations = meanUsedIterations + mover.getUsedIterations()
        maxUsedIterations = max(maxUsedIterations, mover.getUsedIterations())
        meanTime = meanTime + delta * 1000
        if(r):
            meanActualPlateau = meanActualPlateau + mover.getPlateauSize()
            maxActualPlateau = max(maxActualPlateau, mover.getPlateauSize())
            
    meanEdits = meanEdits/len(seeds)
    meanUsedIterations = meanUsedIterations/len(seeds)
    meanTime = meanTime/len(seeds)
    meanActualPlateau = meanActualPlateau/len(seeds)
    return [graph_name, G.numberOfNodes(), getInitName(init), m, s, r, p, 
            meanEdits, minEdits,
            meanUsedIterations, maxUsedIterations,
            meanActualPlateau, maxActualPlateau,
            meanTime]
            
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
                            if(p <= r):
                                df.loc[i] = executeMover(G, name, init, m, s, r, p, b_queue)
                                i += 1
                    else:
                        df.loc[i] = executeMover(G, name, init, m, s, r, 0, b_queue)
                        i+=1
    return df

def writeCSV(file_name, df):
    df['maxIterations'] = df['maxIterations'].apply(np.int64)
    df['plateauSize'] = df['plateauSize'].apply(np.int64)
    df['edits_min'] = df['edits_min'].apply(np.int64)
    df['usedIterations_max'] = df['usedIterations_max'].apply(np.int64)
    df['actualPlateau_max'] = df['actualPlateau_max'].apply(np.int64)
    df['n'] = df['n'].apply(np.int64)
    result_path = output_path + file_name
    df.to_csv(result_path, sep=',', encoding='utf-8')
    
    
    


if not os.path.exists(output_path):
    os.makedirs(output_path)

if graph_set == "small_graphs":
    graphs = ['adjnoun.graph',
          'dolphins.graph', 
          'email.graph', 
          'karate.graph', 
          'lesmis.graph',
          'football.graph',
          'terrorist.edgelist', 
          'grass_web.pairs']
    
if graph_set == "social_network":
    graphs = ['facebook100/Caltech36.mat', 
          'facebook100/Penn94.mat', 
          'lj.edgelist', 
          'orkut.edgelist',
          'youtube.edgelist', 
          'dblp.edgelist', 
          'amazon.edgelist']

if graph_set == "web":
    graphs = ['uk-2002.graph', 
          'eu-2005.graph', 
          'in-2004.graph', 
          'cnr-2000.graph']
    
if graph_set == "biological":
    min_edits = 400
    input_df = pd.read_csv('../input/biological/bio_data.csv')
    graphs = []
    for index, row in input_df.iterrows():
        if(int(row['k']) >= min_edits):
            graphs.append("biological/graphs/" + row['Graph'].split('.')[0] + ".graph")
        #if(len(graphs) == 20):
            #break

    
if graph_set == "generated":
    print('generated')


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
                              'edits_mean',
                              'edits_min',
                              'usedIterations_mean',
                              'usedIterations_max', 
                              'actualPlateau_mean',
                              'actualPlateau_max',
                              'time_mean'])
for graph_name in graphs:
    graph_name_simple = graph_name.split('/')[-1].split('.')[0]
    file_name = graph_name_simple + '_' + scenario +'.csv'
    if(overwrite and os.path.isfile(output_path + file_name)):
        print("Skipping " + graph_name)
        continue
    print("Running QTM on " + graph_name + "...")
    df = runOnGraph(graph_name, df)
    writeCSV(file_name, df)
    print("Finished QTM on " + graph_name)


