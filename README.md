R package for AgeReading
Reference: Storvik, Engebretsen, Aldrin, Fuglebakk: Effect of reader and age reading error on stock assessment: the case of scale readings of Norwegian Spring-Spawning Herring (Clupea harengus)
Folders
  data: Data files
  Mdata: Make R-objects of data used for the routines in exec
  man: Help pages for some routines
  R: General routines
  stan: Stan files
  exec: Files producing results from paper
  
Real data:
  Run_Batch0/1/2.R: Run results on real data on different models
  Produce_results.R: Produce plots.
  
Sensitivity analysis:
  Run_Batch0_sensitivity.R: Run full model on real data, changed hyperparameters
  Run_Batch0_sensitivity.R: Produce plots
  
Simplified simulation study:
  Run_Sim_Batch0.R: Produce results for figures 5 and 6.
  
Simulation study, using posteriors for xi0/xi1:
  Simulation_study_M[1/3/6/Mbar].R: Simulation and fitting using full model 
  Simulation_study_M[1/3/6/Mbar]_indM.R: Simulation and fitting using common M
  Produce_results_sim.R: Produce plots
  
Simulation study 2, also simulating xi0/xi1:
  Simulation2_study_fullM.R: Simulation and fitting using full model 
  Simulation2_study_M1_indM. Simulation and fitting using one M 
  
Figures:
  Produce_results_S2_sim.R:
    - Sim_fit0_xi1est_2.pdf:
  
Files in exec:
  Run_Batch[0/1/2].R: Main analysis on real data, for the three models
  Run_Batch0_sensitivity.R: Main analysis on real data, only full model, sensitivity analysis
  Run_SimRBatch?.R: Simulation 

  
