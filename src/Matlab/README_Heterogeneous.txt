% scripts to test MPL performance on heterogeneous sampling 
%
% Copy Sim_Data, Sim_Analysis, and Sim_Results folders from Zonodo and 
% edit Set_dirNames_MPL_SimData1.txt to provide correct path to Sim_Data, 
% Sim_Analysis, and Sim_Results folders (see Set_dirNames_MPL_SimData1.txt
% for format details).
%
% 1. To generate supplementary Figure 12 from pre-analyzed data
%---------------------------------------------------------------
% Run the following MATLAB scripts in the given order
%
% 1. PreProcessingStep_0_SimData_unEvenTStoch_Loop
% 2. PreProcessingStep_0_SimDataTStoch_Loop
% 3. Plot_Supp_Fig_12
%
% 2. To run complete Heterogeneous analysis simulation: 
%---------------------------------------------------------------
% Run the following MATLAB scripts in the given order
%
% 1. PreProcessingStep_0_SimData_unEvenTStoch_Loop
% 2. PreProcessingStep_0_SimDataTStoch_Loop
% 3. PreProcessingStep_1_SimDataTStoch_Loop
% 4. PreProcessingStep_1_SimData_unEvenTStoch_Loop
% 5. AnalysisMPL_SimDataTStoch_Loop
% 6. AnalysisMPL_SimData_unEvenTStoch_Loop
% 7. Plot_SimData_unEvenTStoch
% 8. Plot_SimDataTStoch
% 9. Plot_SimDataTStoch_diffSelcSites
% 10. Plot_Supp_Fig_12