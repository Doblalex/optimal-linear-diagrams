from datetime import datetime
from collections import defaultdict
import numpy as np
import subprocess
import os
import argparse
import csv
from itertools import chain

def readLD(filepath, invert):
    with open(filepath, 'r', encoding='utf_8') as file:
        csvreader = csv.reader(file)

        header = next(csvreader)

        # file does not have labels for rows and columns
        add_dummy_labels = header[0].strip() in ["0", "1"]

        columns = {}
        rows = {}
        matrix = []

        if not add_dummy_labels:
            for i, column in enumerate(header[1:]):
                columns[i] = column
        else:
            columns = {i: i for i in range(len(header))}

        if not add_dummy_labels:
            i = 0
            for row in csvreader:
                rows[i] = row[0]

                parsedRow = []
                for column in row[1:]:
                    parsedRow.append(int(column))

                matrix.append(parsedRow)

                i += 1
        else:
            i = 0
            for row in chain([header], csvreader):
                rows[i] = i

                parsedRow = []
                for column in row:
                    parsedRow.append(int(column))

                matrix.append(parsedRow)
                i += 1

        matrix = np.matrix(matrix)

        if invert:
            temp = columns
            columns = rows
            rows = temp

            matrix = matrix.T
    return columns, rows, matrix, add_dummy_labels

def count_blocks(matrix: np.ndarray) -> int:
    ans = 0
    for c1, c2 in zip(matrix.T, list(matrix.T[1:])):
        ans += np.count_nonzero(c1 != c2)
    ans += np.count_nonzero(matrix.T[-1])
    ans += np.count_nonzero(matrix.T[0])
    return int(ans/2)

def optimize_gaps_concorde(columns, rows, matrix, force_together=None, performance_tracker=None, row_weights=None):
    starttime = datetime.now()
    if force_together == None:
        force_together = []
    for row in force_together:
        assert(row in rows.values())

    vals_to_columnnames = defaultdict(list)
    MT = matrix.T

    force_together_index = []
    for i, row in rows.items():
        if row in force_together:
            force_together_index.append(i)

    for i, x in columns.items():
        key = tuple(MT[i].tolist()[0])
        vals_to_columnnames[key].append(x)

    columnsets = list()
    index_to_key = dict()
    for i, key in enumerate(vals_to_columnnames):
        columnsets.append(
            set(filter(lambda x: key[x] == 1, list(range(len(key))))))
        index_to_key[i] = key
    # add empty column
    if not any(s == set() for s in columnsets):
        columnsets.append(set())
        index_to_key[len(columnsets)-1] = tuple(0 for _ in range(len(rows)))
    # concorde cannot deal with too small matrices
    while len(columnsets) < 4:
        columnsets.append(set())
        index_to_key[len(columnsets)-1] = None
    force_columns_together = []
    for index in force_together_index:
        columns_together = []
        for i, key in index_to_key.items():
            if key[index] == 1:
                columns_together.append(i)
        force_columns_together.append(columns_together)

    if row_weights != None:
        distance_matrix = []
        for c1 in columnsets:
            distance_matrix.append([])
            for c2 in columnsets:
                dist = sum(row_weights[i] for i in c1.symmetric_difference(c2))
                distance_matrix[-1].append(dist)
                
    else:
        distance_matrix = []
        for c1 in columnsets:
            distance_matrix.append([])
            for c2 in columnsets:
                distance_matrix[-1].append(len(c1.symmetric_difference(c2)))

    count1 = np.count_nonzero(matrix)
    
    for columns_together in force_columns_together:     
        setcolumns = set(columns_together)   
        for i in range(len(distance_matrix)):
            for j in range(len(distance_matrix)):
                if (i in setcolumns) and (j not in setcolumns):
                    distance_matrix[i][j] += (2*count1+1)
                    distance_matrix[j][i] += (2*count1+1)

    strs = []
    strs.append("NAME: TEST")
    strs.append("TYPE: TSP")
    strs.append("DIMENSION: {}".format(len(distance_matrix)))
    strs.append("EDGE_WEIGHT_TYPE: EXPLICIT")
    strs.append("EDGE_WEIGHT_FORMAT: FULL_MATRIX")
    strs.append("EDGE_WEIGHT_SECTION")
    for row in distance_matrix:
        strs.append(' '.join(map(str, row)))
    strs.append("EOF")
    concorde_in = '\n'.join(strs)
    with open("temp.tsp", "w") as f:
        f.write(concorde_in)

    try:
        out = subprocess.check_output(
            ["concorde", "-x", "-o", "solution.sol",
                "temp.tsp"],
            stderr=subprocess.STDOUT).decode("utf-8").split("\n")
    except Exception as e:
        out = e.output.decode("utf-8").split("\n")
    sol = dict()
    for line in out:
        if line.startswith("Optimal Solution:"):
            sol["val"] = float(line.replace(
                "Optimal Solution: ", ""))
        if line.startswith("Total Running Time:"):
            sol["time"] = line.replace("Total Running Time: ", "")
    assert "val" in sol
    with open("solution.sol", "r") as f:
        sol["tour"] = f.read().replace(
            "\n", " ").replace("  ", " ").strip()
        sol["tour"] = list(map(int, sol["tour"].split(" ")))[1:]

    os.system("rm -f solution.sol")
    os.system("rm -f temp.tsp")
    os.system("rm -f temp.sol")

    for i in range(len(sol["tour"])):
        if index_to_key[sol["tour"][i]] == tuple(0 for _ in range(len(rows))):
            start = i
    dist = 0
    for i in range(len(sol["tour"])):
        col1 = columnsets[sol["tour"][i]]
        col2 = columnsets[sol["tour"][(i+1) % len(sol["tour"])]]
        dist += len(col1.symmetric_difference(col2))

    out_M = []
    out_columns = []
    for i in range(start, start+len(sol["tour"])):
        col1 = columnsets[sol["tour"][i % len(sol["tour"])]]
        col2 = columnsets[sol["tour"][(i+1) % len(sol["tour"])]]
        dist += len(col1.symmetric_difference(col2))

        col = index_to_key[sol["tour"][i % len(sol["tour"])]]
        for columnname in vals_to_columnnames[col]:
            out_M.append(list(col))
            out_columns.append(columnname)

    out_columns = {i: x for i, x in enumerate(out_columns)}
    time_taken = datetime.now() - starttime
    if performance_tracker != None:
        performance_tracker[-1].extend(
            [count_blocks(np.array(out_M).T), time_taken.total_seconds()])

    return out_columns, rows, np.array(out_M).T

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", "--input", required=True, type = str)
    parser.add_argument("-out", "--output", required=False, type = str)
    
    return parser.parse_args()
    

def main():
    args = parse_args()
    columns, rows, matrix, dummy_labels = readLD(args.input, False)
    print("Number of blocks before optimization: {}".format(count_blocks(matrix)))
    columns, rows, matrix = optimize_gaps_concorde(columns, rows, matrix)
    print("Number of blocks ater optimization: {}".format(count_blocks(matrix)))
    # example with force_together
    # columns, rows, matrix = optimize_gaps_concorde(columns, rows, matrix, force_together=[rows[0], rows[1]])
    # example with row_weights
    # columns, rows, matrix = optimize_gaps_concorde(columns, rows, matrix, row_weights={i: 1 for i in range(len(rows))})
    if args.output != None:
        with open(args.output, "w") as f:
            writer = csv.writer(f)
            if not dummy_labels:
                writer.writerow([""]+list(columns.values()))
            for i,row in enumerate(rows):
                if not dummy_labels:
                    writer.writerow([rows[row]]+list(matrix[i,:]))
                else:
                    writer.writerow(list(matrix[i,:]))

if __name__ == "__main__":
    main()