# Projective Imaging Toy Model paper
Repository for the code used in the paper Projective Imaging of High-Energy Nuclei via Coherent Exclusive Vector Meson Production in Electron-Nucleus Collisions found here: https://inspirehep.net/literature/2893280
Contains all header files containing the functions to perform the toy model analysis. The macro `FormFactor_Plots_For_Paper.cpp` contain the macros to reproduce the plot in the paper.

To clonde this repo, in terminal: 

```git clone https://github.com/macink/Projective_Imaging_Toy_Model_paper.git```

**Note:** you must have ROOT installed to run the program. Instructions for installation can be found here: https://root.cern/install/

After cloning the repository and installing root, run:

```root -b```

This will start ROOT. You then compile the plotting macro:

```.L FormFactor_Plots_For_Paper.cpp```

**Note:** the header files must be in the same directory that the plotting macro is in for it to compile.
Once compiled, you can produce any of the plots, for example, to produce the first plot, run:

```result_plot()```

This will produce the plot in the same directory that you are currently in.

The macro `FormFactor_tests.cpp` contains an abundance of different plots for analysis regarding the toy model that can be used for reference.
