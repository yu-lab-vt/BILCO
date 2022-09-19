%% Vision Demo

%% Load the problem
load gargoyle-smal.mat

%% Running the minimum-cut algorithm
tic;
[value,cut] = hpf(sim_mat,source,sink);
toc;
% The capacity of the cut is 29696707