# Computing Optimal Linear Diagrams

This repository contains the code for computing optimal linear diagrams using the [Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) TSP solver.

## Requirements

- numpy python package
- The directory containing the [Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) binary has to be added to PATH

## Usage

`python 3 compute_optimal_lineardiagram.py -in INPUT [-out OUTPUT]` where INPUT is the input file containing the linear diagram as csv (either as only 0s and 1s or with column and row labels in the first row and column, respectively) and OUTPUT is the output file name. 

We assume that sets are given as rows and overlaps are given as columns in the input file. This can be changed with the `invert`-parameter of the function `readLD`.

## Forcing sets together

Using the parameter `force_together` of the function `optimize_gaps_concorde`, a list of sets can be forced to be drawn as single line segments, e.g.:

`# columns, rows, matrix = optimize_gaps_concorde(columns, rows, matrix, force_together=[rows[0], rows[1]])`

Note that `rows[0]` and `rows[1]` are the row labels, if no labels were given in the input file, then those are dummy labels.

## Row weights 
Using the parameter `row_weights` of the function `optimize_gaps_concorde` weights can be given to sets in the form of a dictionary where keys are indices of rows and values are integer weights, e.g. (each set has weight 1):

`columns, rows, matrix = optimize_gaps_concorde(columns, rows, matrix, row_weights={i: 1 for i in range(len(rows))})`