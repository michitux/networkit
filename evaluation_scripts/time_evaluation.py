import sys, getopt
import os
import networkit as nk
import timeit
import pandas as pd
import numpy as np
import argparse

#python3 time_evaluation.py -i 3 -m 100 -s True -r True -p 5 -b True

init = 3
maxIterations = 100
sortPaths = True
randomness = True
plateau = 5
bucketQueue_states = [True, False]


graph_sets = ['small_graphs', 'biological400', 'facebook', 'generated', 'social_network', 'web']

def timeMeasurement(graph_set, bucketQueue):
    if graph_set == 'small_graphs':
        graph_names = ['adjnoun.graph', 'dolphins.graph', 'karate.graph', 'lesmis.graph', 'football.graph', 'terrorist.edgelist', 'grassweb.pairs', 'jazz.graph', 'polbooks.graph']
    if graph_set == 'social_network':
        graph_names = ['amazon.edgelist', 'facebook100/Caltech36.mat', 'dblp.edgelist', 'youtube.edgelist', 'facebook100/Penn94.mat', 'lj.edgelist', 'orkut.edgelist']
    if graph_set == 'web':
        graph_names = ['cnr-2000.graph', 'in-2004.graph', 'eu-2005.graph', 'uk-2002.graph']
    if graph_set == 'facebook':
        graph_names = open('facebooknames.txt').read().splitlines()
    if graph_set == 'generated':
        graph_names = open('generatednames.txt').read().splitlines()
    if graph_set == 'biological400':
        graph_names = open('bio400.txt').read().splitlines()

    input_path = '../input/'
    df = pd.DataFrame(columns  = ['graph', 'iteration', 'edits', 'time'])
    row = 0
    for graph_name in graph_names:
        name = graph_name.split('/')[-1].split('.')[0]
        print(name)
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
            df.loc[row] = [name, i, edits[i], time[i]]
            row += 1

    df['iteration'] = df['iteration'].apply(np.int64)
    df['edits'] = df['edits'].apply(np.int64)
    df['time'] = df['time'].apply(np.int64)

    output_path = '../output/QTM_3/' + graph_set + '/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    if bucketQueue:
        file_name = output_path + 'time_bq.csv'
    else:
        file_name = output_path + 'time_lq.csv'

    df.to_csv(file_name, sep=',', encoding='utf-8')


for graph_set in graph_sets:
    for bucketQueue in bucketQueue_states:
        timeMeasurement(graph_set, bucketQueue)
