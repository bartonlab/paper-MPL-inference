%qIn = q(1:100,:);
qIn = q;
lowTh = 1/Nin/2;
highTh = 1 - lowTh;


% definition of trajectory
allSitesAllTrajs = [];
for l = 1:Lin
   thisSiteTraj = qIn(:,l);
   temp1 = thisSiteTraj >= lowTh & thisSiteTraj <= highTh;

   temp = [0 temp1' 0]; % this makes it easy to detect trajectories that have non-zero value at TP =1 and TP = end
   test2 = temp(2:end) - temp(1:end-1);
   thisSiteAllTrajStartInd = find(test2 == 1); 
   thisSiteAllTrajStopInd = find(test2 == -1);
   
   allSitesAllTrajs_Temp = zeros(length(thisSiteAllTrajStartInd),8);
   tCount = 1;
   for tj = 1:length(thisSiteAllTrajStartInd)
       startInd_temp = thisSiteAllTrajStartInd(tj) - 1;% removed -1 for illinworth - 1 ;
       stopInd_temp = thisSiteAllTrajStopInd(tj)-1 + 1;% removed +1 for Illingworth + 1;
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
           allSitesAllTrajs_Temp(tCount,:) = [l startInd_temp stopInd_temp freqAtStartInd_temp freqAtStopInd_temp maxFreq_temp minFreq_temp abs(maxFreq_temp-minFreq_temp)];
           tCount = tCount + 1;
       end
   end
   
   allSitesAllTrajs = [allSitesAllTrajs; allSitesAllTrajs_Temp(1:tCount-1,:)];
      
end

%%


newAllSitesAllTrajs = [];
for l = 1:Lin
    if(sum(allSitesAllTrajs(:,1) == l) > 0)
       allSitesAllTraj_thisSite =  allSitesAllTrajs(allSitesAllTrajs(:,1) == l, :);
       thisSiteTraj = q(:,l);
       newTrajCount = 1;
       newAllSitesAllTrajs_Temp = [];

       currentTrajStartPoint = allSitesAllTraj_thisSite(1,2);
       currentTrajEndPoint = allSitesAllTraj_thisSite(1,3);
       tempTrajThisSite = size(allSitesAllTraj_thisSite,1);
       for tj = 1:tempTrajThisSite

           thisEntryEndPoint = allSitesAllTraj_thisSite(tj, 3);

           % if-else needed to take care of the last entry in allSitesAllTraj_thisSite
           if(tj < tempTrajThisSite)
               nextEntryStartPoint = allSitesAllTraj_thisSite(tj+1, 2);
               condition1 = nextEntryStartPoint - thisEntryEndPoint >= 1;
           else
               nextEntryStartPoint = allSitesAllTraj_thisSite(tj, 2);
               condition1 = true;
           end

           if(condition1) 
               % this means more than 2 TP difference between the end and
               % start of the two trajs 

               % thisTraj has ended: update currentTrajEndPoint and check
               % if currentTrajStartit meets condition to retain trajcetory
               currentTrajEndPoint = thisEntryEndPoint;


               maxFreq_temp = max(thisSiteTraj(currentTrajStartPoint:currentTrajEndPoint));
               minFreq_temp = min(thisSiteTraj(currentTrajStartPoint:currentTrajEndPoint));

               % make entry of new traj ONLY if it is greater that noise
               % threshold of thisSiteThresh1
               thisSiteThresh1 = 0.10;
               if(abs(maxFreq_temp-minFreq_temp) > thisSiteThresh1)               
                   newAllSitesAllTrajs_Temp(newTrajCount,:) = [l currentTrajStartPoint currentTrajEndPoint thisSiteTraj(currentTrajStartPoint) thisSiteTraj(currentTrajEndPoint) maxFreq_temp minFreq_temp abs(maxFreq_temp-minFreq_temp)];
                   newTrajCount = newTrajCount + 1;
               end

               if(tj < tempTrajThisSite)
                   currentTrajStartPoint = allSitesAllTraj_thisSite(tj+1,2);
                   currentTrajEndPoint = allSitesAllTraj_thisSite(tj+1,3);
               end
           else
               % join prevTraj and thisTraj
               % do nothing, joining will be done in the next iter of tj
           end
       end

       newAllSitesAllTrajs = [newAllSitesAllTrajs; newAllSitesAllTrajs_Temp];
    end
end

siteTrajAll = zeros(500,2,300);
trajPerSite = zeros(1,Lin);
for l = 1:Lin
    indTrajThisSites_logical = newAllSitesAllTrajs(:,1) == l;
    trajPerSite(l) = sum(indTrajThisSites_logical);
    siteTrajAll(1:trajPerSite(l),:,l) = newAllSitesAllTrajs(indTrajThisSites_logical,2:3);
end


%%

% first make a temp qTemp1 matrix which only replaces 'sampling noice' 0's
% and 1's with 1/N and 1 - 1/N respectively
qTemp1 = q;
for l = 1:Lin
    newAllSitesAllTrajs_thisSite = newAllSitesAllTrajs(newAllSitesAllTrajs(:,1) == l, :);
    for tj = 1:size(newAllSitesAllTrajs_thisSite,1)
        thisEntryStartPointPlus1 = newAllSitesAllTrajs_thisSite(tj,2) + 1;
        thisEntryEndPointMinus1 = newAllSitesAllTrajs_thisSite(tj,3) - 1;
        
        qThisPartial = qTemp1(thisEntryStartPointPlus1:thisEntryEndPointMinus1,l);
        qThisPartial(qThisPartial == 0) = 1/Nin;
        qThisPartial(qThisPartial == 1) = 1 - 1/Nin;
        qTemp1(thisEntryStartPointPlus1:thisEntryEndPointMinus1,l) = qThisPartial;
    end
end


% now make new q trajectories by selection only the trajectory time points
% in newAllSitesAllTrajs

% this loop only updats the *non trajectory* points
qTemp3 = zeros(size(q,1), size(q,2));
for l = 1:Lin
    newAllSitesAllTrajs_thisSite = newAllSitesAllTrajs(newAllSitesAllTrajs(:,1) == l, :);
    if(size(newAllSitesAllTrajs_thisSite,1) > 0)
        for tj = 1:size(newAllSitesAllTrajs_thisSite,1)
            thisEntryStartPointPlus1 = newAllSitesAllTrajs_thisSite(tj,2);
            thisEntryEndPointMinus1 = newAllSitesAllTrajs_thisSite(tj,3);
            % make firstSegment 
            if(tj == 1)
                firstSegmentStart = 1;
            else
                firstSegmentStart = newAllSitesAllTrajs_thisSite(tj-1,3);
            end
            firstSegmentEnd = newAllSitesAllTrajs_thisSite(tj,2);
            firstSemantValue = round(qTemp1(firstSegmentEnd, l));

            qTemp3(firstSegmentStart:firstSegmentEnd,l) = firstSemantValue;

            if(tj == size(newAllSitesAllTrajs_thisSite,1))
                secondSegmentStart = newAllSitesAllTrajs_thisSite(tj,3);
                secondSegmentEnd = numSamplingPoints;
                secondSegmantValue = round(qTemp1(secondSegmentEnd, l));
                qTemp3(secondSegmentStart:secondSegmentEnd,l) = secondSegmantValue;
            end
        end
    else
        % for sites that do not have any entry in newAllSitesAllTrajs_thisSite
        temp = mean(q(:,l)); % this value will be ither close to 1 or 0 depending upon the trajectory
        qTemp3(:,l) = round(temp);
    end
end

%%
% this loop only updats the *trajectory* points
for l = 1:Lin
    newAllSitesAllTrajs_thisSite = newAllSitesAllTrajs(newAllSitesAllTrajs(:,1) == l, :);
    for tj = 1:size(newAllSitesAllTrajs_thisSite,1)
        thisEntryStartPoint = newAllSitesAllTrajs_thisSite(tj,2);
        thisEntryEndPoint = newAllSitesAllTrajs_thisSite(tj,3);
        qTemp3(thisEntryStartPoint:thisEntryEndPoint,l) = qTemp1(thisEntryStartPoint:thisEntryEndPoint,l);
    end
end

validQ = double(qTemp3 ~= 0 & qTemp3 ~= 1);
validSelCoeff = validQ;

