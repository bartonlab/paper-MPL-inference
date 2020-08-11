%%

function synOrNonSynBeforeFilt = markSynNonSyn(numNT, numSitesNT, numSitesAA, numUniqueTimePoints, codonCompListCell, countOfCodonsPerSiteAA, msaREFCodonCell)
%  0: syn
%  1: non-syn
% -1: can not be classified as syn/non-syn, no variation at this site to decide 
%  2: poly site, where at least 1 NT is a gap

% clasifies only those NTs that are observed in the DATA
synOrNonSynBeforeFilt = -10*ones(1, numNT*numSitesNT);
for i = 1:numSitesAA

    REFCodonAllTPThisAASite = zeros(numUniqueTimePoints, 3);
    allCodonsThisSiteTemp = [];
    tpVecTemp = [];
    for tp = 1:numUniqueTimePoints
        allCodonsThisSiteTemp = [allCodonsThisSiteTemp; codonCompListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i),:)];
        tpVecTemp = [tpVecTemp; tp*ones(size(codonCompListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i),:), 1),1)];
        REFCodonAllTPThisAASite(tp,:) = msaREFCodonCell{tp,i};
    end
    
    
    allCodonsThisSite = unique(allCodonsThisSiteTemp, 'rows');
    for k = 1:3 % loop over 3 NT indicis (positions) of the codon
        uniqueNTThisSite = unique(allCodonsThisSite(:,k));
        numPolyNTatKthPosition = length(uniqueNTThisSite);
        if(numPolyNTatKthPosition > 1)
            for l = 1:numNT
                thisNTSiteNonSyn = -1;
                thisNTInt = l;
                if(l == 5)
                    thsiNTInt = 16;
                end
                
                % do only if l is a memebr of uniqueNTThisSite, check if it
                % is syn or nonSyn
                if(ismember(l, uniqueNTThisSite))
                    thisNTSiteNonSyn = 0;
                    for tp = 1:numUniqueTimePoints
                        REFCodon = REFCodonAllTPThisAASite(tp,:);
                        REFAA = nt2aa(int2nt(REFCodon));
                        codonsThisSiteThisTP = allCodonsThisSiteTemp(tpVecTemp == tp,:);
%                         REFCodon
%                         codonsThisSiteThisTP
%                         pause
                        numCodonsThisSiteThisTP = size(codonsThisSiteThisTP, 1);
                        % run tests only if this TP has more than 1 codons
                        if(numCodonsThisSiteThisTP > 1)
                            if( l ~= REFCodon(k)) 
                                % select all codons from the list
                                % codonsThisSiteThisTP that have NT 'l' at
                                % position k  
                                temp30 = codonsThisSiteThisTP(:,k) == l; 

                                selectedCodons = codonsThisSiteThisTP(temp30,:);
                                %[i k l tp]
                                %pause
                                for yy = 1:sum(temp30)
%                                     if(i == 30)
%                                         [i k l tp]
%                                         nt2aa(int2nt(selectedCodons(yy,:)))
%                                         [REFAA nt2aa(int2nt(selectedCodons(yy,:)))]
                                         temp40 = unique([REFAA nt2aa(int2nt(selectedCodons(yy,:)), 'AlternativeStartCodons', false)]);
%                                         pause
%                                     end
                                    numberOFAAs = length(temp40);
                                    if(numberOFAAs > 1 && thisNTSiteNonSyn ~= 1)
                                        thisNTSiteNonSyn = 1;
                                    end
                                end
                            end
                        end
                        
                    end
                end
                
                synOrNonSynBeforeFilt((i-1)*3*numNT + (k-1)*numNT + l) = thisNTSiteNonSyn;
            end
        else
            synOrNonSynBeforeFilt((i-1)*3*numNT + (k-1)*numNT + (1:numNT)) = -1;
        end
    end
end