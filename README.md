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

main: set up a world made of a voxel list a vasculature. Add healthy stroma cells in each voxel and add a few clonogenic tumor cells in the center. Creates a healthy vasculature.
Then run the class simulaton that take for argument a list of process and the predefine world.



## Shell Script to run the simulation

This script takes three input arguments: the name of the Python input file ($INFILE), the number of iterations to run ($ITER), and the name of the configuration file ($CONFIG_NAME). It automatically determines the starting count number based on the existing directories in the ./output/$CONFIG_NAME directory. It then creates a new directory for each iteration of the simulation and copies the input files ($INFILE and all Python files) and support files to the new directory. It also generates a new random seed for each iteration and writes it to the input file.
For each iteration, the script generates a new job submission script ($SCRIPT) using the cat command and submits it to the LSF job scheduler using the bsub command. The job name is set to "$INFILE-$CONFIG_NAME-$COUNT", where $INFILE is the name of the Python input file, $CONFIG_NAME is the name of the configuration file, and $COUNT is the iteration number. The script then increments the count and repeats the process until the specified number of iterations is reached.

./run_AMBER.sh main.py 10 my_config

## 


