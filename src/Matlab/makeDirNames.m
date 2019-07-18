function [] = makeDirNames()
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    disp('Error: system si not unix and not PC...')
    pause
end
currDir = pwd;
dirNameData = [currDir(1:end-6) 'wfsim' chosenSlash 'data' chosenSlash];
dirNameAnalysis = [currDir chosenSlash 'Analysis' chosenSlash];
if(ispc)
    dirNameABC = [currDir(1:end-6) 'external' chosenSlash 'WFABC_v1.1' chosenSlash 'binaries'  chosenSlash 'Windows'  chosenSlash];
elseif(isunix)
    dirNameABC = [currDir(1:end-6) 'external' chosenSlash 'WFABC_v1.1' chosenSlash 'binaries'  chosenSlash 'Linux'  chosenSlash];
end
dirNameApproxWF = [currDir(1:end-6) 'external' chosenSlash 'Approxwf' chosenSlash];

fileID = fopen('dirNames.txt','w');
fprintf(fileID,'%s\n', '% This file specifies directory paths to Data, Analysis, ABC, ApproxWF');
fprintf(fileID,'%s\n', '% the syntax is /dir1/dir2/dir3/ ');
fprintf(fileID,'%s\n', '% for windows based system, the code will automatically reformat the path');
fprintf(fileID,'%s\n', '%');
fprintf(fileID,'%s\n', '% Please fill in the ''/.../'' with the actual path of your installation.');
fprintf(fileID,'%s\n', '%');
fprintf(fileID,'%s\n', '% Note: dirNameABC for Linux based system /.../paper-MPL-inference-master/src/external/WFABC_v1.1/binaries/Linux/');
fprintf(fileID,'%s\n', '%       dirNameABC for Windows based system C:\...\paper-MPL-inference-master\src\external\WFABC_v1.1\binaries\Windows\');
fprintf(fileID,'%s\n', '% -------------------------------------------------------------------------');
fprintf(fileID,'%s\n', ['dirNameData=' dirNameData]);
fprintf(fileID,'%s\n', ['dirNameAnalysis=' dirNameAnalysis]);
fprintf(fileID,'%s\n', ['dirNameABC=' dirNameABC]);
fprintf(fileID,'%s\n', ['dirNameApproxWF=' dirNameApproxWF]);
fclose(fileID);

