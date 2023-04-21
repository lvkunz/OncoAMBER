# AMBER (Agent-based Modeling of Biophysics Evolution for Radiotherapy)
## Introduction

In our approach, we use a hybrid model that combines agent-based components with some continuum-like aspects.
This allows us to create efficient multiscale simulations using Monte Carlo techniques.
Our model uses tumor properties and radiation effects on individual cells as input to predict macro-scale tumor evolution.
We work with a voxelized geometry where cells divide based on a representative statistical distribution and diffuse to neighboring voxels depending on the crowding in each voxel.
To provide oxygen to the cells, we utilize a pre-generated 'healthy' vasculature.
We've also conducted sub-voxel size simulations to determine a representative distribution of a cell's oxygen within a voxel, depending on crowding and vessel density.
This helps to mimic chronic hypoxia. As the tumor grows, it can block and destroy the vasculature, leading to acute hypoxia.
In response, hypoxia triggers VEGF production and subsequent angiogenesis towards the tumor.
The angiogenesis process is based on a directed random walk that takes crowding and the VEGF gradient into account.

High crowding or pressure can lead to long-lasting hypoxia of cells and eventually result in their necrosis. 
We'll integrate the radiation effects on individual cells using the simple interface provided by TOPAS and TOPAS-nBio, as well as other in-house DNA damage and DNA repair models.

## Structure

import amber 

amber.World(--)
amber.Cell(--)
amber.Voxel(--)
amber.Vessel(--)
amber.VesselNetwork(--)
amber.Simulator(--)
etc. 

example on github
