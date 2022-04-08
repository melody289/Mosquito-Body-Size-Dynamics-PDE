# Mosquito-Body-Size-Dynamics-PDE

The following files run various codes 
dataFormat.m 
pulls the data from file (access this https://data.mendeley.com/datasets/2ddbwns4gb/1) and formats it into a matrix count

difCombVar.m 
Will run the 13 parameter situations for simulating both high and low density, This file creates the data set labeled ParCombFit.mat.

difHighVar.m 
Will run the 13 parameter situations for simulating only high density data, This file creates the data set labeled ParHighFit.mat.

difLowVar.m 
Will run the 13 parameter situations for simulating only low density data, This file creates the data set labeled ParLowFit.mat.

ParError.m
This will give the error for all simulations run. This needs the data file from https://data.mendeley.com/datasets/2ddbwns4gb/1 and it will need a parameter Anewg and/or Anewgl.
ch= 0 is low density run only, load ParLowFit.mat.
ch= 1 is high density run only, load ParHighFit.mat.
ch= 2 is both low and high density run, load ParCombFit.mat.


C_DecayMass.m
This runs the simulations with decaying resources which is included in the growth rate of the mosquito larvae.


C_DecayDeath.m
This runs the simulations with decaying resources which is included in the death rate of the mosquito larvae.


C_CyclicMass.m
This runs the simulations with cyclic resources which is included in the growth rate of the mosquito larvae.

C_CyclicDeath.m
This runs the simulations with cyclic resources which is included in the death rate of the mosquito larvae.

CyclicMassCombine.m
This was to take all the files saved from C_CyclicMass.m into a single matrix.

Note: that several file runs were too large to include.


IC2.mat Contains inital conditions for some of the runs.

The following is the code already prerun.
VarSepFit.mat  # For Figure 3
VarCombFit.mat  # For Figure 4
cyclicgroup.mat # For Figure 9 and 10
