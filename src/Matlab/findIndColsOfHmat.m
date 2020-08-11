
function [selcSitesEpi, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(Hmat)

rrefMat = rref(Hmat);

rrefMat = round(1000*rrefMat)/1000;
k1 = 1;
k2 = 1;
clustCount = 1;
condRepeat = 1;
potential = [];
cluster = [];
while(condRepeat == 1)
   % if the (k1,k2) element is a 1, then this column aly
   % corresponds to a site that can be estimated uniquely
   if(rrefMat(k1,k2) == 1)
       % the corresponding (column) site can be estimated uniquely IF the
       % rest of the row is all zeros
       if(sum(rrefMat(k1,k2+1:end) ~= 0) == 0)
           potential = [potential k2];
       % but if the rest of the rows contains any non zero entry, the non
       % zero entry column(s) and k2th column form a cluster. We can only
       % estimate the 'sum' of SC of this clustr, not individual SC.
       else
           tempCluster = [k2 find(rrefMat(k1,k2+1:end) ~= 0) + k2];
           % IF this is the 1st cluster, just add it to the cluster list
           if(clustCount == 1)
               cluster{clustCount} = tempCluster;
               clustCount = clustCount + 1;
           % IF this is NOT the 1st cluster, check if the tempCluster has
           % any overlap with any of the existing cluster. If they have
           % overlap, join it with that cluster, otherwie make a new
           % cluster entry
           else
               i = 1;
               clusterOverlapFound = 0; 
               % check if this cluster has ovrlap with any of the existing
               % clusters
               while(i <= clustCount-1 && clusterOverlapFound == 0)
                   interTemp = intersect(cluster{i},tempCluster);
                   if(~isempty(interTemp))
                       clusterOverlapFound = 1;
                   end
                   i = i + 1;
               end
               
               % update existing cluser
               if(clusterOverlapFound == 1)
                   tempNew = union(cluster{i-1}, tempCluster);
                   cluster{i-1} = tempNew;
               % add new cluster
               else
                   cluster{clustCount} = tempCluster;
                   clustCount = clustCount + 1;
               end
           end
           
       end
       k1 = k1 + 1;
       k2 = k2 + 1;
   % if (k1,k2) element is zero, move to next column in the same row    
   elseif(rrefMat(k1,k2) == 0)
       k2 = k2 + 1;
   else
       disp(['rrefMat(k1,k2) = ' num2str(rrefMat(k1,k2)) '...'])
       pause
   end
   
   if(k1 > size(rrefMat,1) || k2 > size(rrefMat,1))
       condRepeat = 0;
   end
end

wellCondSites = potential;
illCondSites = setdiff(1:size(rrefMat,1), wellCondSites);

selcSitesEpi = logical(zeros(1, size(rrefMat,1)));
selcSitesEpi(wellCondSites) = true;