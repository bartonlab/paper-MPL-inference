%% Code to filter trajecotries for ABC method

% allSitesAllTrajs contain 
% [site startIndOfTraj stopIndOfTraj freqAtStartIndOfTraj freqAtStopIndAtTraj maxFreqInThisTraf minFreqInThisTraj]

qIn = q;
lowTh = 1/Nin/2;
highTh = 1 - lowTh;


% definition of trajectory
allSitesAllTrajs = [];
for l = 1:L
   thisSiteTraj = qIn(:,l);
   temp1 = thisSiteTraj >= lowTh & thisSiteTraj <= highTh;

   temp = [0 temp1' 0]; % this makes it easy to detect trajectories that have non-zero value at TP =1 and TP = end
   test2 = temp(2:end) - temp(1:end-1);
   thisSiteAllTrajStartInd = find(test2 == 1); 
   thisSiteAllTrajStopInd = find(test2 == -1);
   
   allSitesAllTrajs_Temp = zeros(length(thisSiteAllTrajStartInd),7);
   tCount = 1;
   for tj = 1:length(thisSiteAllTrajStartInd)
       startInd_temp = thisSiteAllTrajStartInd(tj) - 1;
       stopInd_temp = thisSiteAllTrajStopInd(tj)-1 + 1;
       if(startInd_temp < 1)
           startInd_temp = 1;
       end
       if(startInd_temp > size(qIn,1))
           startInd_temp = size(qIn,1);
       end
       
       if(stopInd_temp < 1)
           stopInd_temp = 1;
       end
       if(stopInd_temp > size(qIn,1))
           stopInd_temp = size(qIn,1);
       end
       
       
       
       freqAtStartInd_temp = thisSiteTraj(startInd_temp);
       freqAtStopInd_temp = thisSiteTraj(stopInd_temp);
       maxFreq_temp = max(thisSiteTraj(startInd_temp:stopInd_temp));
       minFreq_temp = min(thisSiteTraj(startInd_temp:stopInd_temp));
       if(stopInd_temp - startInd_temp + 1 > 0)
           allSitesAllTrajs_Temp(tCount,:) = [l startInd_temp stopInd_temp freqAtStartInd_temp freqAtStopInd_temp maxFreq_temp minFreq_temp];
           tCount = tCount + 1;
       end
   end
   
   allSitesAllTrajs = [allSitesAllTrajs; allSitesAllTrajs_Temp(1:tCount-1,:)];
      
end

