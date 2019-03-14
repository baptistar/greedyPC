# Polynomial Chaos Randomized Greedy Algorithm (RGA)

## What is the RG Algorithm?

The RGA is a greedy approach for building a sparse Polynomial Chaos (PC) approximation of a function given input and output sample evaluations. At each iteration, the algorithm evaluates a random subset of basis functions from a large candidate dictionary to greedily select terms to add to the approximation without having a computational cost that scales with the dictionary size. This is combined with an efficient QR-based strategy to update the expansion and a leave-one-out estimate to measure generalization error. For more information on the details of the algorithm, please see the [JCP paper](https://www.sciencedirect.com/science/article/pii/S0021999119300865).

## Authors

Ricardo Baptista (MIT) and Prasanth Nair (University of Toronto)

E-mails: <rsb@mit.edu> or <ricarsb@gmail.com>, <pbn@utias.utoronto.ca>

## Installation

The Randomized Greedy Algorithm is implemented in MATLAB. This package includes all necessary classes and methods to run the RGA. To compare to PC approximations that are constructed using L1 minimization, the code requires the [SPGl1 library](https://www.cs.ubc.ca/~mpf/spgl1/) to be available in the local path. To assess the accuracy of the PC approximations, the scripts evaluate the polynomial models at Sparse Grid and pseudo-random Sobol points. To generate these points, the code requires the [Sparse Grid](https://people.sc.fsu.edu/~jburkardt/m_src/sparse_grid_cc/sparse_grid_cc.html) and [Sobol](https://people.sc.fsu.edu/~jburkardt/m_src/sobol/sobol.html) libraries to be available in the local path.

## Example running RGA on a benchmark problem

We provide an example of running the RGA on an algebraic test problem with 5 dimensional inputs. The PC approximation is built using training data ranging in size from N = 10 to N = 1000 input points and function evaluations. The code compares the RGA to a complete basis PC expansion, L1 minimization, and two variants of Orthogonal Matching Pursuit. The code can also be run from MATLAB using the file `scripts/example_run.m`.

The script first defines the input parameters. These include the dimension of the inputs `d`, the basis functions `basis`, the maximum order of the total degree polynomial basis `order`, and the anonymous target function `func`.

For each N, we sample a uniform set of points on ![equation](https://latex.codecogs.com/gif.latex?%5B-1%2C1%5D%5E%7Bd%7D) and evaluate the function at these inputs. Only these sample points are seen by the algorithm during the training phase.
	
	Xall = 2*rand(N,d) - 1;
	Yall = func(Xall);

For each sample size, N, the function extracts a subset of the training points `X,Y` from `Xall,Yall` and builds the sparse PC model using `fit`. Given a set of testing points `XTest, YTest, wTest`, the `TestErr` function estimates the relative error in the mean and standard deviation as well as the relative mean-squared error of the approximation.

	RGA = RandomizedGreedy(d, order, basis);
	RGA = RGA.fit(X,Y);
	[TestErr, MeanErr, StdDevErr] = RGA.TestErr(XTest, YTest, wTest);

The script produces a plot with the convergence results from running RGA and other sparse PC algorithms for increasing N.
