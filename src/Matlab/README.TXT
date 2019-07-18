%==========================================================================
%                                 README.TXT
%==========================================================================
%
%
% This file contains information on running the provided scripts
% 
%
% Dependencies: MATLAB R 2017a 9.2.0 or higher 
%
% MATLAB toolbox required: Global Optimization Toolbox 3.4.2 or higher
%   
%
%
% Running WFABC
%-----------------------------------------
%  1. Download the WFABCv1.1 package from http://jjensenlab.org/software
%  2. Unpack the zip file in folder src/external/
%     Verify the directory structure after installation is
%     src/external/WFABC_v1.1/binaries/
%  3. Run script DataGenABC_comp.m to covnert ground truth (GT) data into 
%     format readable by WFABC method.
%  4. Run script AnalysisABC_comp.m to run the WFABC analysis (Foll, 2015) 
%     on GT data.
%  5. Step 4 will generate .csv file containing estimates of selection 
%    coefficients and save the file in the directory src/Matlab/
%
%
% Running ApproxWF
%-----------------------------------------
%  1. Download the ApproxWF package (commit fcc7964 dated: 2016-10-04) from 
%     https://bitbucket.org/wegmannlab/approxwf/wiki/Home
%  2. Unack the zip file in folder src/external/Approxwf
%     Verify the directory structure after installation is
%     src/external/Approxwf
%     NOTE: Verify the Version of Approxwf is fcc7964. This script may not 
%     work properly with later versions of the Approxwf package.
%  3. Run script DataGenApproxWF_comp.m to covnert GT data into format
%     readable by ApproxWF method.
%  4. Run script AnalysisApproxWF_comp.m to run the ApproxWF analysis 
%     (Feder, 2014) on GT data.
%  5. Step 4 will generate .csv file containing estimates of selection 
%     coefficients and save the file in the directory src/Matlab/
%
%
% Running FIT
%-----------------------------------------
%  1. Run script AnalysisFIT_comp.m to run the FIT analysis (Feder, 2014) 
%     on GT data.
%  2. Step 1 will generate .csv file containing estimates of selection 
%     coefficients and save the file in the directory src/Matlab/
%
%
% Running IM
%-----------------------------------------
%  1. Run script AnalysisIM_comp.m to run the method of Illingowrth
%     (Illingworth, 2011) on GT data.
%  2. Step 1 will generate .csv file containing estimates of selection 
%     coefficients and save the file in the directory src/Matlab/
%
%  Note: The script AnalysisIM_comp.m requires a platform dependent MEX 
%        function (func1_v3_mex.mexa64 / func1_v3_mex.mexw64). The current 
%        version of the code works for
%            (a)  64-bit Linux system
%            (b)  64-bit Windows system
%