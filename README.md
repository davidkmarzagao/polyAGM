# polyAGM
Package for Property Prediction of Polymers based on Automatically Generated Motifs. We use graph kernel together with probabilistic methods to infer property predictions and their explanations.

This work was developed by Shannon R. Petersen, David Kohan Marzagão, Georgina L. Gregory, Yichen Huang, David A.
Clifton, Charlotte K. Williams, Clive R. Siviour at the departments of Engineering, Chemistry and Computer Science of the University of Oxford. 

The experiments described in the paper (and fully reproducible here) were run on a M1 Pro MacBook Pro, 32GB. The requirements can be installed with 'pip install -r requirements_polyagm.txt' (or just run "installing_packages.ipynb"). The time to install and run experiments should be minimal. If you are evaluating feature vectors using rdkit, each polymer may take a few miuntes to compute. 

To repruduce results shown on Figure 5 of the manuscript, just execute cells in order in file "[Figure 5] - Parity Plots.ipynb". You may want to change parameters of function run_kernels_against_target in order to get extra results. For step-by-step procedure, see "predictions.ipynb". For a demo on how SMILES are tranformed into graphs, please see notebook "From SMILES to Graphs.ipynb". 