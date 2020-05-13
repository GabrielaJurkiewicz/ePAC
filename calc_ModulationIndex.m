function [MI,MISurr,MeanAmp] = calc_ModulationIndex(Map,MessMap,MeanLowFreq,nbBins,PhaseBins,Nboot,p,plotWithMask)

    MI = zeros(size(Map,1),1);
    MISurr = zeros(Nboot,size(Map,1));
    MeanAmp = zeros(size(Map,1),nbBins); 
    dPhi    = 2*pi/nbBins;
    
    if sum(sum(Map~=0))>0

        %% assign amplitudes to Phase Bins
        Phase   = angle(hilbert(MeanLowFreq));
        MeanAmp = zeros(size(Map,1),nbBins); 
        MeanAmpSurr = zeros(Nboot,size(Map,1),nbBins); 
        for b = 1:nbBins   
            idx = find((Phase <  PhaseBins(b)+dPhi) & (Phase >=  PhaseBins(b)));
            MeanAmp(:,b) = mean(Map(:,idx),2); 
            
            if plotWithMask
                for nb = 1:Nboot
                    MeanAmpSurr(nb,:,b) = squeeze(mean(MessMap(nb,:,idx),3));
                end
            end
        end
        
        MeanAmp = MeanAmp./repmat(sum(MeanAmp,2),1,size(MeanAmp,2));
        MeanAmp(isnan(MeanAmp)) = 0;
        MI = (log(nbBins)-(-sum(MeanAmp.*log(MeanAmp),2)))/log(nbBins);
        MI(isnan(MI)) = 0;

        if plotWithMask
            
            MeanAmpSurr = MeanAmpSurr./repmat(sum(MeanAmpSurr,3),1,1,size(MeanAmpSurr,3));
            MeanAmpSurr(isnan(MeanAmpSurr)) = 0;
            MISurr = (log(nbBins)-(-sum(MeanAmpSurr.*log(MeanAmpSurr),3)))/log(nbBins);
            MISurr(isnan(MISurr)) = 0;
            
            %% threshold MeanAmp values using Exteme Stat
            for fa = 1:size(Map,1)
                thresh = prctile(squeeze(max(MeanAmpSurr(:,fa,:),[],3)),p);
                tmp = MeanAmp(fa,:);
                tmp(tmp<thresh) = 0;
                MeanAmp(fa,:) = tmp;
            end
            
            %% Normalize MI using mean of MISurr
            tmp = squeeze(mean(MISurr,1));
            MI = MI-tmp';
            for nb = 1:Nboot
                MISurr(nb,:) = squeeze(MISurr(nb,:))-tmp;
            end
        
        end
    end
end