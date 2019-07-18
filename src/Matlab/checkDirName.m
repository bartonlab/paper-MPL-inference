%% this code formats dir/filename.ext for use in linux command prompt. Linux
% requires that if a directory name has a space, it should have a front
% slash in serted in it
%                                |          |     |
%                                v          v     v
% /local/staff/ee/mssohail/Matlab\ Codes/HIV\ Time\ Series/P3/F1/SAM_files/reads_p3_1_F1_sorted.sam
%
%
%
% Last updated Oct 2016
% Author: M Saqib Sohail

% ---------------------------------------
function correctedSourceSAMDirFileName = checkDirName(sourceSAMDirFileName)

correctedSourceSAMDirFileName = sourceSAMDirFileName;
indOfSpace = find(sourceSAMDirFileName == ' ');
for i = length(indOfSpace):-1:1
    correctedSourceSAMDirFileName = [correctedSourceSAMDirFileName(1:indOfSpace(i)-1) '\' correctedSourceSAMDirFileName(indOfSpace(i):end)];
end
