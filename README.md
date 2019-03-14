### Polynomial Chaos Randomized Greedy Algorithm (RGA)

# What is the RG Algorithm?

The RGA is a greedy approach for building a sparse Polynomial Chaos approximation of a function given input and output sample points. At each iteration, the algorithm evaluates a random subset of basis functions from a large candidate dictionary to greedily select terms to add to the approximation without having having a computational cost that scales with the dictionary size. This is combined an efficient QR-based strategy to update the expansion and a leave-one-out error estimate to measure generalization error. For more information on the details of the algorithm, please see the JCP paper.

# Authors

Ricardo Baptista (MIT) and Prasanth Nair (University of Toronto)
Emails: rsb@mit.edu or ricarsb@gmail.com, pbn@utias.utoronto.ca

# Installation

The Greedy PC algorithm is implemented in MATLAB. The code includes all necessary packages to run the RGA,. In order to assess the accuracy of the model, we generate a set of using a Sparse Grid package. For the test problems, we also compare to using Quasi-Monte Carlo by evaluating the points at a set of .  We also compare to the L1 regularization using the SPGL1 library.

-> i4_sobol package
-> SPGL1 in the path
-> sparse grid package 

# Example running RGA on a benchmark problem

We provide an example for running the RGA on an algebraic test problem with 5 dimensionsal inputs. The PC approximation is built with input samples and function evaluates of size N = 10 to N = 1000. The code compares the RGA to a Full PC expansion (FPC), L1 Minimization (BPDN), and two variants of Orthogonal Matching Pursuit (OMP). The code can also be run from MATLAB using the file scripts/example_run.m 

The script first defines the input parameters. These include the dimension of the inputs, the basis of the Polynomials.

We randomnly sample a uniform set of points from the function as well

For each sample size, N, the function extracts a subset of the training points and builds the PC model 

RGA = RandomizedGreedy(d, order, basis);
RGA = RGA.fit(X,Y);
[TestE, MeanE, StdE] = RGA.TestErr(xTest, Ytest, wTest);

The results (mean, standard deviation and test error) are plotted.
