![Title](AMBER_title.png)

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

## Installation

Download the github repository and run the following command in the main directory:

```bash
bash install.sh
```

This will install the latest version of the amber package itself. Make sure to activate your python environment before running the command.

github: https://github.com/lvkunz/OncoAMBER; PyPi: https://pypi.org/project/OncoAMBER/

## Usage

I recommand specifying the path to the amber package at the beginning of your python script:

example:
```python
import sys
sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
```

you can then import the amber package.

```python
import amber

amber.World(-)
amber.Cell(-)
amber.Voxel(-)
amber.Vessel(-)
amber.VesselNetwork(-)
amber.Simulator(-)
etc
```

## Running the example locally

The example.py file is explained in further details in the Tutorial.ipnyb notebook.
However, you need to know that amber needs a CONFIG file to be run. This file contains all the parameters needed to run the simulation.

you need to use the following terminal command to run the example:
```bash
python example.py CONFIG
```

Note that the CONFIG does not include the .txt extension! 

If you're running the example through a GUI editor such as PyCharm, you can specify the CONFIG file in the run configuration. 

On PyCharm, in the top right corner, click on the dropdown menu and select "Edit Configurations...".
In the "Script parameters" field, type the name of the CONFIG file you want to use.

## Running the example on the cluster

To run the example on the cluster, you need to use the following command:

```bash
bash run_AMBER.sh example.py <n_iter> CONFIG
```
You might need to change the shell script for your own cluster.

Where <n_iter> is the number of iterations you want to run. You can then decide which parameters you want to change in the CONFIG file for each iteration.
For example when prompted 
```
Enter a parameter name to change, or type 'none' to continue: 
```

you can reply 'dt' to change the time step for each iteration. Choose all the parameters you want to change, then type 'none' to continue.
Then you will be prompted to enter the values for each parameter you chose to change for each iteration.

The script will then run all simulations and store them in an output folder. You'll find that the folder hierarchy goes through dates, name of config file and pyhton
script, and then the iteration number. The output folder will contain all the data from the simulation, including the CONFIG file used for each iteration.
The output folder will also contain a log file with the date and time of the simulation, as well as the parameters used for each iteration stored in a .csv file.
If you'd like to change the seed to some random value I recommand setting it to -1 in the CONFIG file. The script will then generate a random seed for each iteration.

To easily read the results of your simulations you can run the plotting,py script by changing the output folder directory you want to use and the name of the parameters you want to look at.
You can also change t_min, t_max and choose to plot the results for all iterations or just a few of them.

