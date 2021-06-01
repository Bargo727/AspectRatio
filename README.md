# AspectRatio
This repository contains all code used to output figures pertaining to the lattice model (LM) in the paper "Using Cell Shape to Control Spatiotemporal Dynamics of Synthetic Microbial Consortia".  The name of each file corresponds to the figure in the paper for which the indicated code was used to produce.  All are coded in MATLAB.

#fig4 contains files used to produce figure 4 in the main paper.  fig4a.m is used to produce the dynamics of the nematic ordering parameter q versus the parameter p_{rot} for various trap sizes.  Follow the comments within the script to obtain figures as seen in the paper.  fig4bc provides code for the time series snapshots of invasion and bulk displacement mechanisms in the lattice model.  Read the commentary in the code.  You can modify the initial conditions of the code by uncommenting appropriate code lines.  The code will also provide code for computing average dynamics.  fig4de.m is used to visualize the dynamics of bulk displacement and invasion.  You can create movies from that code.  Also, following the commentary will allow you to produce fig 4d in the main paper.  You can use the same file to produce the supplemental figures that pertain to lattice model bullk force and displacement.  Commented code will also let you produce the figure in the supplement that shows probabilities of horizontal growth for each strain in an AB wall format. 

#suppfig_master_equation.m produces figures that pertain to the master equation seen in the supplementary information.  Modify code within according to comments to obtain the necessary figures.

#videos contains videos that illustrate bulk displacement, invasion, and spatial oscillations in the LM.  The name of the video files describes what the video file depicts.  This folder contains three files:  Spatial_Oscillation.mp4, Bulk_Displacment.mp4, and Invasion.mp4
