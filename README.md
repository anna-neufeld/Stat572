# Stat572
All code for Stat 572 paper on "Low Risk Estimators of Population Size". 

- The files atom_mcmc.R and beta_mcmc.R implement Johndrow et al.'s Bayesian fitting procedures for the case where H is a 3-atom discrete mixture and the case when H is Beta(a,b), respectively.

- Run_Violin_Plots.R calls atom_mcmc and beta_mcmc and conducts the main replication of Johndrow et al.'s paper. It creates the results plotted in the 10-panel violin plot for simulated data and the 2-panel violin plot for the snowshoe hare data. 
    - Results from running this file are stored in myOutput_513 folder
    - To actually create the plots using the stored results, use Make_Violon_Plots.R
    - Raw PNGs of the plots are also stored in myOutput_513

- Run_Optimal_Alpha_Simulation.R has all the code for the experiment where I plot bias, variance, and MSE as a function of alpha for 3 different scenarios. In this file, the M_h model is fit with Maximum Likelihood, not a Bayesian algorithm. 



- Analyze_Extra_PointEstimate_Sims.R stores
