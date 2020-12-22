# AspectRatio
This repository contains all code used to output figures pertaining to the lattice model in the paper "Using Cell Shape to Control Spatiotemporal Dynamics of Synthetic Microbial Consortia".  The name of each file corresponds to the figure in the paper for which the indicated code was used to produce.  All are coded in MATLAB.

#fig5 contains files used to produce figure 5 in the main paper.  fig5a.m is used to produce the dynamics of the nematic ordering parameter q versus the parameter p_{rot} for various trap sizes.  Follow the comments within the script to obtain figures as seen in the paper.  fig5bc provides code for the time series snapshots of invasion and bulk displacement mechanisms in the lattice model.  Read the commentary in the code.  You can modify the initial conditions of the code by uncommenting appropriate code lines.  The code will also provide code for computing average dynamics.  fig5de.m is used to visualize the dynamics of bulk displacement and invasion.  You can create movies from that code.  Also, following the commentary will allow you to produce fig 5d in the main paper.  

#fig6 contains files to produce figure 6 in the main paper.  fig6ae provides code to see spatial oscillations in the lattice model; within it is also commented code that allows one to compute the average strain fraction dynamics as is seen in fig 6e.  Follow the code commentary.  fig 6bcd provides code that produces oscillations from the reduced model in the paper, and, if you follow the commentary, will show you code for producing the plots the show the dependence of frequency on t_d and amplitude on H.  The remaining file, EffectiveOscillation.m, is a script that has the equations for the reduced model.

#suppfig_master_equation.m produces figures that pertain to the master equation seen in the supplementary information.  Modify code within according to comments to obtain the necessary figures.
