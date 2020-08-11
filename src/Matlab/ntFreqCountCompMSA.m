function freqCount = ntFreqCountCompMSA(msa_int, seq_freq)

% Similar to ntFreqCount but modified for compressed MSA where the MSA only
% contains unique sequences and a separate vector of frequeices for each
% uniue sequence

% Returns the frequency of each amino acid at each protien position
% input: MSA in integer format --  (sequences) x (protien length)
% output: freqCount ((17) x (protien length) frequency count of each amino
% acid at each protien location. 

% use nt2int to convert MSA to integer values. Each column of the returnd 
% freqCount matix indicates the frequency count of NT at that position
% in the protien. The frequency count of 'A' is given at the index 
% corresponding to its NT integer value.

% 17th entry corresponds ot the code 0 which is for (*)

% find frequency of amino acids at each position in the protien 
freqCount = zeros(17, size(msa_int, 2)); % 17 as MATLAB uses 17 int to represent NT


for i = 1:size(msa_int,2)
    for currentNT = 1:16
        thisNTvec = msa_int(:,i) == currentNT; % contains a 1 where currentNT exists in this vector, 0 otherwise
        freqCount(currentNT, i) = sum(seq_freq(thisNTvec));
    end
    thisNTvec = msa_int(:,i) == 0 | msa_int(:,i) == 17; % contains a 1 where currentNT exists in this vector, 0 otherwise
    freqCount(17, i) = sum(seq_freq(thisNTvec));
    
end

return