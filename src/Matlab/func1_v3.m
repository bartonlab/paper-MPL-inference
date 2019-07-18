%% calculate sigmaEst_i_eff

% this is v3, removes the making of siteTrajAll from this function and
% places it in the outer code...its now passed inside as a 3D array


% v2 modifiedwhich modifiies the num, den and adds temp90 for the forard
% and back ward equation


%function [overAllLogLikelihood, logLikelihood]= func1(sigmaEst ,numSamplingPoints, L, timePointsToDrop, lastValidTimePoint, qijAll, validQ, validSelCoeff, dT, q, ng)

function overAllLogLikelihood = func1_v3(sigmaEst, numSamplingPoints, L, timePointsToDrop, lastValidTimePoint, dT, ng, trajPerSite, validQ, validSelCoeff, q, siteTrajAll, gMat)
%function overAllLogLikelihood = func1_v3(sigmaEst ,numSamplingPoints, L, timePointsToDrop, lastValidTimePoint, validQ, validSelCoeff, dT, q, ng, siteTrajAll, trajPerSite, gMat)
%#codegen

%fprintf('calculating effective selection coefficent (sigmaEff)...')
sigmaEst_eff = zeros((numSamplingPoints - timePointsToDrop), L);

for t = 1:(numSamplingPoints - timePointsToDrop)
%     for i = 1:L
%         sigmaEst_eff(t, i) = sigmaEst(i) + sigmaEst*gMat(i,:,t)';
%     end
    sigmaEst_eff(t,:) = sigmaEst' + gMat(:,:,t)*sigmaEst';
    
end

% for t = 1:(numSamplingPoints - timePointsToDrop)
% %    t
%      qijAtTimeTk = qijAll{t};
%      for i = 1:L
%          % if qi < threshold, sigmaEst_eff = 0, pp. 992
%          if(validQ(t,i) == 1) 
%         %display('y')
%              tempSum = 0;
%              for j = 1:L
%                  
%                 % qj < threshold, ignore its linkage effect pp. 992
%                 if(~(i==j) && (validSelCoeff(t,j) == 1)) 
%                     thisQ = qijAtTimeTk{i,j};
%                     
%                     q00 = thisQ(1,1);
%                     q01 = thisQ(1,2);
%                     q10 = thisQ(2,1);
%                     q11 = thisQ(2,2);
%                     
%                     
%                     if(q00 == 0 && q01 == 0)
%                         display(' Error 0: q00 and q01 are both zero!')
%                         [q00 q01; q10 q11]
%                         i
%                         j
%                         t
%                         pause
%                     end
%                     if(q11 == 0 && q10 == 0)
%                         display(' Error 1: q11 and q10 are both zero!')
%                         [q00 q01; q10 q11]
%                         i
%                         j
%                         t
%                         pause
%                     end
%                     
% %                     if(i == 2)
% %                         tempSum
% %                     end
%                     tempSum = tempSum + sigmaEst(j)*(q11/(q11 + q10) - q01/(q01 + q00));
% %                     if(i == 2)
% %                             [i j t]
% %                             [q00 q01; q10 q11]
% %                             sigmaEst(j)*(q11/(q11 + q10) - q01/(q01 + q00))
% %                             tempSum
% %                             pause
% %                     end
%                         %display('YY')
%                 end
%              end
%              sigmaEst_eff(t, i) = sigmaEst(i) + tempSum;
%          end
%      end
% end
%display('done')


%% %% plot the deterministic trajectory based on eq 8, 10, and estimate of
% sigma_i 



% display('----------------------------')
%             display('curent sigma...')
%              sigmaEst
%             pause

deltaTk = dT;%10;
qEst = 0*ones(numSamplingPoints, L);
logLikelihood = 0*ones(numSamplingPoints, L);
for l = 1:L
    thisSite = l;
    temp1 = siteTrajAll(:,:,thisSite);
    if(~isempty(temp1)) % do if this site has atleast 1 trajectory
        numTrajThisSite = trajPerSite(l);% size(temp1,1);
        for k = 1:numTrajThisSite
            thisTrajStartEnd = temp1(k,:);
            trajStartPoint = thisTrajStartEnd(1);
            trajEndPoint = thisTrajStartEnd(2);
            
            trajMiddlePointQc = trajStartPoint + ceil((trajEndPoint - trajStartPoint)/2);
            qEst(trajMiddlePointQc,thisSite) = q(trajMiddlePointQc,thisSite);

            for t = (trajMiddlePointQc + 1):trajEndPoint
               temp90 = exp(sigmaEst_eff(t-1,l)*deltaTk);
               if(isinf(temp90))
%                    display('in inf...')
                  if(temp90 > 0)
                      num = exp(700); 
                  else
                      num = exp(-700); 
                  end
               else
                  num = qEst(t-1,l)*temp90;
               end
               
               if(isinf(temp90))
                   if(temp90 > 0)
                       den = exp(700);
                   else
                       den = exp(-700);
                   end
               else
                  den = (1 - qEst(t-1,l) + qEst(t-1,l)*temp90); 
               end
               
               
               qEst(t,l) = num/den;
               % equation 11
               this_qObs = q(t,thisSite)*validQ(t,thisSite);
               
               
               firstTerm = nchoosek(ng, round(ng*this_qObs));
               likelihood = firstTerm * ( (qEst(t,l))^(ng*this_qObs) ) * ( (1-qEst(t,l))^(ng*(1-this_qObs)) ); 
               
               logLikelihood(t,l) = log10(likelihood);
               if(isinf(log10(likelihood)))
                   logLikelihood(t,l) = -320;
               end
               if(isnan(logLikelihood(t,l)))
                   fprintf('Forward...')
                   [l t num den num/den q(t,l) likelihood logLikelihood(t,l)]
                   sigmaEst
                   pause
               end
               %[l t num den num/den q(t,l) likelihood logLikelihood(t,l)]
               
            end
            
            for t = trajMiddlePointQc -1:-1:trajStartPoint
               num = qEst(t+1,l);
               
               temp92 = exp(sigmaEst_eff(t+1,l)*deltaTk);
               if(isinf(temp92))
%                    display('in inf...')
                   if(temp92 > 0)
                      den = exp(700); 
                   else
                       den = exp(-700); 
                   end
               else
                  den = (temp92 - qEst(t+1,l)*temp92 + qEst(t+1,l));
               end
               
               %num/den
               qEst(t,l) = num/den;
               %pause
               
               % equation 11
               this_qObs = q(t,thisSite)*validQ(t,thisSite);
               firstTerm = nchoosek(ng, round(ng*this_qObs));
               likelihood = firstTerm * ( (qEst(t,l))^(ng*this_qObs) ) * ( (1-qEst(t,l))^(ng*(1-this_qObs)) ); 
               logLikelihood(t,l) = log10(likelihood);
               if(isinf(log10(likelihood)))
                   logLikelihood(t,l) = -320;
               end
               %[l t num den num/den q(t,l) likelihood logLikelihood(t,l)]
               if(isnan(logLikelihood(t,l)))
                   fprintf('Backward...')
                   [l t num den num/den q(t,l) likelihood logLikelihood(t,l)]
                   pause
               end
            end
        end        
    end
%     figure
%     plot(qEst(1:lastValidTimePoint,thisSite), 'k.')
%     hold on
%     plot(q(1:lastValidTimePoint,thisSite).*validQ(1:lastValidTimePoint,thisSite))%(qOrig(:,thisSite))
end

overAllLogLikelihood = -sum(sum(logLikelihood));

%sigmaEst
%overAllLogLikelihood
%display('----------------------------------------------------------------')
% %%
% 
% % plot tranjectory at each site
% for l = 1:4
%     figure(l)
%    subplot(2,1,1)
%    plot(validQ(1:lastValidTimePoint,l).*q(1:lastValidTimePoint,l))
% %    axis([3000 3500 -1 1])
%    xlabel('time points')
%    ylabel('Frequency Count')
%    subplot(2,1,2)
%    plot(validSelCoeff(1:lastValidTimePoint,l).*sigmaEst_eff(1:lastValidTimePoint,l))
%    %plot(logLikelihood(1:lastValidTimePoint,l))
%    xlabel('time points')
%    ylabel('Effective selection')
%    title(['Site number: ' num2str(l)])
% % axis([3000 3500 -10*10^-3 10*10^-3])   
% end
