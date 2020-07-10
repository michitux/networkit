import sys, getopt
import os
import networkit as nk
import datetime
import pandas as pd
import numpy as np


global graph_set
global scenario

try:
   opts, args = getopt.getopt(sys.argv[1:],"g:s:",["graph_set=","scenario="])
except getopt.GetoptError:
   print ('evaluation.py -g <graph_set> -s <scenario> ')
   sys.exit(2)

for opt, arg in opts:
    if opt in ("-g", "--graph_set"):
        graph_set = arg
    elif opt in ("-s", "--scenario"):
        scenario = arg
        
global seeds
seeds = [i for i in range(9)]
global input_path
input_path = "../input/"
global output_path 
output_path = '../output/QTM/' + graph_set + '/'
        
        
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
        a = datetime.datetime.now()  
        mover.run()   
        delta = datetime.datetime.now() - a
        meanEdits = meanEdits + mover.getNumberOfEdits()
        minEdits = min(minEdits, mover.getNumberOfEdits())
        meanUsedIterations = meanUsedIterations + mover.getUsedIterations()
        maxUsedIterations = max(maxUsedIterations, mover.getUsedIterations())
        meanTime = meanTime + delta.total_seconds() * 1000
        if(r):
            meanActualPlateau = meanActualPlateau + mover.getPlateauSize()
            maxActualPlateau = max(maxActualPlateau, mover.getPlateauSize())
            
    seed_count = len(seeds)*1.0
    meanEdits = meanEdits/seed_count
    meanUsedIterations = meanUsedIterations/seed_count
    meanTime = meanTime/seed_count
    meanActualPlateau = meanActualPlateau/seed_count
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
        G = nk.readGraph(graph_path, nk.Format.EdgeListTabOne)
    if(graph_name.split('.')[1] == "pairs"):
        G = nk.readGraph(graph_path, nk.Format.EdgeListTabOne)
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

def writeCSV(file_name, df):
    df['maxIterations'] = df['maxIterations'].apply(np.int64)
    df['plateauSize'] = df['plateauSize'].apply(np.int64)
    df['edits_min'] = df['edits_min'].apply(np.int64)
    df['usedIterations_max'] = df['usedIterations_max'].apply(np.int64)
    df['actualPlateau_max'] = df['actualPlateau_max'].apply(np.int64)
    df['n'] = df['n'].apply(np.int64)
    df = df.round({'edits_mean':2, 'usedIterations_mean': 2, 'actualPlateau_mean': 2, 'time_mean':3})
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

global Intializations
global maxIterations
global sortPaths
global randomness
global plateauSize
global b_queue


if(scenario == 'full'):
    maxIterations = [0, 5, 100]
    initializations = [0, 1, 2, 3]
    sortPaths = [True, False]
    randomness = [True, False]
    plateauSize = [2, 5, 10]
    b_queue = True
if(scenario == 'plateauBound'):
    initializations = [0, 1, 2, 3]
    maxIterations = [100]
    sortPaths = [True, False]
    randomness = [True]
    plateauSize = [100]
    b_queue = True
if(scenario == 'withoutBucketQueue'):
    initializations = [0, 1, 2, 3]
    maxIterations = [0, 5, 100]
    sortPaths = [True, False]
    randomness = [True, False]
    plateauSize = [2, 5, 10]
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
file_name = graph_set + '_' + scenario +'.csv'
for graph_name in graphs:
    print(graph_name)
    df = runOnGraph(graph_name, df)
writeCSV(file_name, df)


