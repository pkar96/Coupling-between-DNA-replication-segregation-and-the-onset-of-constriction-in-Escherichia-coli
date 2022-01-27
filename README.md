# Coupling-between-DNA-replication-segregation-and-the-onset-of-constriction-in-Escherichia-coli
Contains MATLAB codes used to generate Figure 3 and Figure 5 of the manuscript, "Coupling between DNA replication, segregation and the onset of constriction in Escherichia coli".

Code used to generate Figure 3 is named as Figure_3.m.

Code used to generate Figure 5 is present in Figure 5 folder. The folder contains files-
1. ao_model_change_c.m - Simulates the adder per origin model
2. ch_model_change_c.m - Simulates the Cooper Helmstetter model
3. parallel_adder_change_c.m - Simulates the parallel adder model
4. indep_adder_change_c.m - Simulates the independent adder model
5. call.m - Calls one of the above cell cycle models. Output obtained contains cell cycle characteristics such as C+D period, generation time, length at birth, division and initiation for multiple cells in a population. These cell cycle characteristics can be plotted against time from thymine shift. 
6. Cell.m - A class which creates an object of type Cell. The object contains all the cell charcteristics required during the cell cycle models' simulation.  
7. binning_with_error_1.m - Used to get the binned data given the raw data and the bin edges.
