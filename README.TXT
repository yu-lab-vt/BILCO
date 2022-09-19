This package implements bidirectional pushing with linear component operations (BILCO) in C++ language and provides a MATLAB wrapper for it.
-----------------------------
PACKAGE
-----------------------------
BILCO
-----------------------------
  ./BILCOMex.cpp - the C++ code of the wrapper

  ./BILCO.m - function to build the wrapper

  ./BILCOMex.mexw64 - binary file for the MEX-function compiled using MATLAB R2020a + MSVC 2019

  ./DTW_Edge_input.cpp - the C++ code of the DTW with inputting distance matrix

  ./DTW_Edge_input.mexw64 - binary file for DTW_Edge_input.cpp compiled using MATLAB R2020a + MSVC 2019

USING BILCO
-----------------------------

[minCut,maxFlow] = BILCO(ref,tst,Gij,smo,initialCut,winSize,metric)

Output:

  minCut - the height of min-cut in each column. Can be converted to the path through midPoint2path(minCut)

  maxFlow - max-flow value

Input: (N: number of time series pairs, T1: length of reference curve, T2: length of test curve)
 
  ref - reference curve, N x T1
 
  tst - test curve, N x T2

  Gij - relationship between time series pairs. For example, [1,2] means 1st pair and 2nd pair are related

  smo - hyperparameter kappa/2, balance the profile similarity and distance between related warping functions

  initialCut - initialization information  N x (T1-1)

  winSize - the allowed window size for misalignment. Default is 10000

  metric - metric for distance measure. 'squared'(default) or 'absolute'

Experiment
-----------------------------
  ./demo_syn.m - compare methods on synthetic data

  ./real_data_application1.m - compare methods on application to calculating signal propagation

  ./real_data_application2.m - compare methods on application to extracting depth information

  ./real_data_application3.m - compare methods on application to signature identification

OTHER PACKAGES
-----------------------------
  ./referenced_methods - implementation of IBFS, HPF, BK, HIPR with MATLAB wrapper
