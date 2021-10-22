# Installation instructions
#### To run the experiments, you need to install the rerf package, which computes the geodesic forests, and locally change one of the files (as described below) to run our extended version of GF.
0) In R run ''install.packages("remotes")''
1) Load rerf_2.0.4 package here: https://cran.r-project.org/src/contrib/Archive/rerf/
2) Go to the downloaded package and replace the file "rerf/R/UrerfHelperFunctions.R" with the corresponding file in the "code" folder.
3) Go to R and run
library("remotes")
install_local("path to rerf folder")

# Run G-KSG
The stand alone G-KSG method is implemented in "code/G_KSG.R" and can be called with its default parameters given X,Y (both matrices) and k (typically equal to 3). For an example on how to apply the method, see "code/run_example.R". To execute the code, it is assumed that you are in the "code" direktory.


# Run Experiments
The files for running experiments are in the folder code/experiments.

To execute the code, you must, however, be in the directory "code" (not in "code/experiments").

The experiments for increasing number of samples, which document the runtime are in file "evaluate_runtime_and_accuracy_per_sample.R" and all other experiments can be executed from "run_test.R" by setting the described parameters accordingly.
