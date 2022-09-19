%% Demo General
% This demo generates the following similrity matrixand 
% finds the corresponding graph minimum-cut
clear;
clc;

%% Generating the similarity matrix
a = sparse([0 3.1 3.2 0 0 0
    0 0 0 2 0 0
    0 0 0 1.5 2 0
    0 0 0 0 0 1.3
    0 0 0 0 0 3.3
    0 0 0 0 0 0 ]);

%% SOURCE node is node# 1
source = 1;

%% SINK node is node# 6
sink = 6;

% a = a*100000000;
%% Running the minimumcut algorithm
tic;
[value,cut] = hpf2(a,source,sink);
toc;
% The capacity of the cut is 5 and Nodes 1 (the source) and 2 are the
% source set

% a = sparse([0,2.1;0,0]);
% [value,cut] = hpf2(a,1,2)