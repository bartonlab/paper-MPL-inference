
function [aucLinkItrTemp, aucLinkSelcItrTemp, ...
    MCitrsUsable_Pos_itr, MCitrsUsable_Neg_itr] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutLink, sigmaEstOutLinkSelc, ...
                                                            posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp)

if(sum(posOnlyItrTemp) == 0 || isempty(posOnlyItrTemp))
    aucLinkItrTemp(1) = -1;
else
    [~,~,~, aucLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink, 1);
end

if(sum(negOnlyItrTemp) == 0 || isempty(negOnlyItrTemp))
    aucLinkItrTemp(2) = -1;
else
    [~,~,~, aucLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink, 1);
end

% --- calculations for filtered ---
if(sum(posOnlySelcItrTemp) == 0 || isempty(posOnlySelcItrTemp))
    aucLinkSelcItrTemp(1) = -1;
elseif(sum(posOnlySelcItrTemp) == length(posOnlySelcItrTemp)) % add a pseudo neutral value to true and predicted
    MCitrsUsable_Pos_itr = 1;%0;
    [~,~,~, aucLinkSelcItrTemp(1)] = perfcurve(double([posOnlySelcItrTemp 0]), [sigmaEstOutLinkSelc; 0], 1);
else
    [~,~,~, aucLinkSelcItrTemp(1)] = perfcurve(double(posOnlySelcItrTemp), sigmaEstOutLinkSelc, 1);
end

if(sum(negOnlySelcItrTemp) == 0 || isempty(negOnlySelcItrTemp))
    aucLinkSelcItrTemp(2) = -1;
elseif(sum(negOnlySelcItrTemp) == length(negOnlySelcItrTemp)) % add a pseudo neutral value to true and predicted
    MCitrsUsable_Neg_itr = 1;
    [~,~,~, aucLinkSelcItrTemp(2)] = perfcurve(double(~logical([negOnlySelcItrTemp 0])), [sigmaEstOutLinkSelc; 0], 1);
else
    [~,~,~, aucLinkSelcItrTemp(2)] = perfcurve(double(~negOnlySelcItrTemp), sigmaEstOutLinkSelc, 1);
end


% 
% function [aucLinkItrTemp, aucLinkSelcItrTemp, ...
%     MCitrsUsable_Pos_itr, MCitrsUsable_Neg_itr] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutLink, sigmaEstOutLinkSelc, ...
%                                                             posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp)
% 
% if(isempty(posOnlyItrTemp))
%     aucLinkItrTemp(1) = -1;
% else
%     [~,~,~, aucLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink, 1);
% end
% 
% if(isempty(negOnlyItrTemp))
%     aucLinkItrTemp(2) = -1;
% else
%     [~,~,~, aucLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink, 1);
% end
% 
% % --- calculations for filtered ---
% 
% % MCitrsUsable_Pos_itr = 1;
% % if(isempty(posOnlySelcItrTemp))
% %     MCitrsUsable_Pos_itr = 1;%0;
% % elseif(sum(posOnlySelcItrTemp) == length(posOnlySelcItrTemp)) % add a pseudo neutral value to true and predicted
% %     MCitrsUsable_Pos_itr = 1;%0;
% %     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 0]);
% %     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 0]);
% %     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; 0];
% % elseif(sum(posOnlySelcItrTemp) == 0) % adda pseudo positive value to true and predicted
% %     MCitrsUsable_Pos_itr = 1;%0;
% %     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 1]);
% %     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 0]);
% %     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; max(perSiteSelction)];
% % end
% 
% 
% % MCitrsUsable_Neg_itr = 1;
% % if(isempty(negOnlySelcItrTemp))
% %     MCitrsUsable_Neg_itr = 1;%0;
% % elseif(sum(negOnlySelcItrTemp) == length(negOnlySelcItrTemp)) % add a pseudo neutral value to true and predicted
% %     MCitrsUsable_Neg_itr = 1;%0;
% %     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 0]);
% %     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 0]);
% %     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; 0];
% % elseif(sum(negOnlySelcItrTemp) == 0) % add a pseudo negative value to true and predicted
% %     MCitrsUsable_Neg_itr = 1;%0;
% %     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 1]);
% %     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 0]);
% %     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; min(perSiteSelction)];
% % end
% 
% 
% if(sum(posOnlySelcItrTemp) == length(posOnlySelcItrTemp)) % add a pseudo neutral value to true and predicted
%     MCitrsUsable_Pos_itr = 1;%0;
%     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 0]);
%     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 0]);
%     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; 0];
% elseif(sum(posOnlySelcItrTemp) == 0 && sum(negOnlySelcItrTemp) ~= 0)%(sum(posOnlySelcItrTemp) == 0) % adda pseudo positive value to true and predicted
%     MCitrsUsable_Pos_itr = 1;%0;
%     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 1]);
%     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 0]);
%     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; max(perSiteSelction)];    
% end
% 
% if(sum(negOnlySelcItrTemp) == length(negOnlySelcItrTemp)) % add a pseudo neutral value to true and predicted
%     MCitrsUsable_Neg_itr = 1;%0;
%     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 0]);
%     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 0]);
%     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; 0];
% elseif(sum(negOnlySelcItrTemp) == 0 && sum(posOnlySelcItrTemp) ~= 0)%(sum(negOnlySelcItrTemp) == 0) % add a pseudo negative value to true and predicted
%     MCitrsUsable_Neg_itr = 1;%0;
%     negOnlySelcItrTemp = logical([negOnlySelcItrTemp 1]);
%     posOnlySelcItrTemp = logical([posOnlySelcItrTemp 0]);
%     sigmaEstOutLinkSelc = [sigmaEstOutLinkSelc; min(perSiteSelction)];
% end
% 
% % if(MCitrsUsable_Pos_itr == 1)
% %     [~,~,~, aucLinkSelcItrTemp(1)] = perfcurve(double(posOnlySelcItrTemp), sigmaEstOutLinkSelc, 1);
% % end
% % if(MCitrsUsable_Neg_itr == 1)
% %     [~,~,~, aucLinkSelcItrTemp(2)] = perfcurve(double(~negOnlySelcItrTemp), sigmaEstOutLinkSelc, 1);
% % end
% 
% 
% 
% if(sum(posOnlySelcItrTemp) == 0)
%     aucLinkSelcItrTemp(1) = -1;
% else
%     [~,~,~, aucLinkSelcItrTemp(1)] = perfcurve(double(posOnlySelcItrTemp), sigmaEstOutLinkSelc, 1);
% end
% 
% if(sum(negOnlySelcItrTemp) == 0)
%     aucLinkSelcItrTemp(2) = -1;
% else
%     [~,~,~, aucLinkSelcItrTemp(2)] = perfcurve(double(~negOnlySelcItrTemp), sigmaEstOutLinkSelc, 1);
% end
% 
