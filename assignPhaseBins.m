function [meanAmp,meanAmpSurr] = assignPhaseBins(Map,MessMap,MeanLowFreq,nbBins,PhaseBins,Nboot,plotWithMask)

    meanAmp = zeros(size(Map,1),nbBins);
    meanAmpSurr = zeros(Nboot,size(Map,1),nbBins);
    dPhi        = 2*pi/nbBins;
    
    if sum(sum(Map~=0))>0

        %% assign amplitudes to Phase Bins
        Phase   = angle(hilbert(MeanLowFreq)); 
        for b = 1:nbBins   
            idx = find((Phase <  PhaseBins(b)+dPhi) & (Phase >=  PhaseBins(b)));
            meanAmp(:,b) = mean(Map(:,idx),2); 
            
            if plotWithMask
                for nb = 1:Nboot
                    meanAmpSurr(nb,:,b) = squeeze(mean(MessMap(nb,:,idx),3));
                end
            end
        end
        
    end
    
end