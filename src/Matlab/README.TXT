%==========================================================================
%                                 README.TXT
%==========================================================================
%
%
% This file contains information on running the provided scripts
% 
%
% Dependencies           : MATLAB R 2017a 9.2.0 or higher 
%
% MATLAB toolbox required: Global Optimization Toolbox 3.4.2 or higher
%   
%
% Preliminaries 
%-----------------------------------------
%  1. install WFABC and ApproxWF from thir respective sources.
%  2. Set paths for data and analysis directories in dirNames.txt
%  3. The script AnalysisIM_comp.m requires a platform dependent MEX 
%     function (func1_v3_mex.mexa64 / func1_v3_mex.mexw64). The current 
%     version of the code works for
%          (a)  64-bit Linux system
%          (b)  64-bit Windows system
%
%
% Scripts and their functions
%-----------------------------------------
%
% 1. Script DataGenABC_comp.m covnerts GT data into format readable by
%     WFABC method
% 2. Script DataGenApproxWF_comp.m covnerts GT data into format
%    readable by ApproxWF method.
% 3. Script AnalysisABC_comp.m runs the WFABC analysis (Foll, 2015).
% 4. Script AnalysisApproxWF_comp.m runs the ApproxWF analysis (Ferrer, 2016). 
% 5. Script AnalysisFIT_comp.m runs the FIT analysis (Feder, 2014).
% 6. Script AnalysisIM_comp.m runs the Deterministic analysis of
%    Illingworth (Illingworth, 2011). 
% 7. Script Plot_comp.m plots a comparison of the above methods and also
%    saves the analysis results in .CSV format
%
% Running scripts for various datasets
%-----------------------------------------
% 
% 1. Set the path in dirName.txt
% 2. Specify the dataset to analyse in each script by setting the
%    appropriate value for the variable 'thisSet'.
% 3. Run the script
%
% Note: Analysis parameters for each dataset are pre-specified in
%       getSysParam_4.m