#!/usr/bin/env python3

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert protein .cm files to metis graph files')
    parser.add_argument('--threshold', help='The threshold', default=0, type=float)
    parser.add_argument('input', help='The input file')
    parser.add_argument('output', help='The output file')

    args = parser.parse_args()

    n = m = 0

    with open(args.input, 'r') as input_file:
        for ln, line in enumerate(input_file):
            if ln == 0:
                n = int(line)
                neighbors = [[] for i in range(n)]
            elif ln <= n:
                continue
            elif line.rstrip():
                u = ln - n
                
                expected_neighbors = n - u

                all_neighbors = line.split('\t')

                assert(len(all_neighbors) == expected_neighbors)

                for i, weight in enumerate(map(float, all_neighbors)):
                    v = u + i + 1

                    if weight >= args.threshold:
                        neighbors[u-1].append(v)
                        neighbors[v-1].append(u)
                        m += 1

    with open(args.output, 'w') as output_file:
        print("{} {} 0".format(n, m), file=output_file)
        for neigh in neighbors:
            print("{}".format(" ".join(map(str, neigh))), file=output_file)
