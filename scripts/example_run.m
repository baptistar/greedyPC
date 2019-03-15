% Polynomial Chaos Randomized Greedy Algorithm
%
% Copyright (C) 2019 R. Baptista & P. Nair
%
% GreedyPC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% GreedyPC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License 
% along with GreedyPC. If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2019 MIT & University of Toronto
% Authors: Ricardo Baptista & Prasanth Nair
% E-mails: rsb@mit.edu or ricarsb@gmail.com & pbn@utias.utoronto.ca
%

% Test script for Least Squares PCE code
clear; close all; clc

addpath('../SpectralToolbox')
addpath('../TestProblems')
addpath('../Methods')
addpath('../tools')

% define test parameters
d = 5;
order = 3;
grid_level = 5;
basis = 'Legendre';
func = @(x) algebraic_1(x);
N = floor(logspace(1,3,10));

% generate training data: uniform inputs on [-1,1]
Xall = 2*rand(max(N), d) - 1;
Yall = func(Xall);

% call sparse grid to determine the test samples
[XTest, wTest] = sparse_grid(d, grid_level, 'unif');
YTest = func(XTest);

% declare vectors to store results
TestE = zeros(length(N),5);
MeanE = zeros(length(N),5);
StdE  = zeros(length(N),5);

for i=1:length(N)
    
    fprintf('Test: N = %d\n', N(i));
    
    % extract samples
    X = Xall(1:N(i),:);
    Y = Yall(1:N(i));
    
    % train PC model using LeastSquares 
    FPC  = LeastSquaresPCE(d, order, basis);
    FPC  = FPC.fit(X,Y);
    [TestE(i,1), MeanE(i,1), StdE(i,1)] = FPC.testErr(XTest, YTest, wTest);

    % train PC model using BPDN 
    BPDN = L1Minimization(d, order, basis);
    BPDN = BPDN.fit(X,Y);
    [TestE(i,2), MeanE(i,2), StdE(i,2)] = BPDN.testErr(XTest, YTest, wTest);
    
    % train PC model using traditional OMP 
    OMP  = OrthogonalMatchingPursuit(d, order, basis, 'old');
    OMP  = OMP.fit(X,Y);
    [TestE(i,3), MeanE(i,3), StdE(i,3)] = OMP.testErr(XTest, YTest, wTest);

    % train PC model using OMP with optimal basis selection strategy
    OMPN = OrthogonalMatchingPursuit(d, order, basis, 'opt');
    OMPN = OMPN.fit(X,Y);
    [TestE(i,4), MeanE(i,4), StdE(i,4)] = OMPN.testErr(XTest, YTest, wTest);

    % train PC model using Randomized Greedy algorithm
    RGA  = RandomizedGreedy(d, order, basis);
    RGA  = RGA.fit(X,Y);
    [TestE(i,5), MeanE(i,5), StdE(i,5)] = RGA.testErr(XTest, YTest, wTest);

end

% plot results
figure('position',[0,0,1200,600]);

subplot(1,3,1)
loglog(N, MeanE, '-o', 'MarkerSize', 10)
xlabel('$N$')
ylabel('Mean Error')
legend({'FPC','BPDN','OMP','OMPN','RGA'})

subplot(1,3,2)
loglog(N, StdE, '-o', 'MarkerSize', 10)
xlabel('$N$')
ylabel('Std Error')
legend({'FPC','BPDN','OMP','OMPN','RGA'})

subplot(1,3,3)
loglog(N, TestE, '-o', 'MarkerSize', 10)
xlabel('$N$')
ylabel('Test Set Error')
legend({'FPC','BPDN','OMP','OMPN','RGA'})

% -- END OF FILE --