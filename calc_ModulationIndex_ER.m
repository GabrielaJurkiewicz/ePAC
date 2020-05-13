function [MI,MXS,MeanAmp] = calc_ModulationIndex_ER(MeanAmp,MeanAmpSurr,nbBins,Nboot,p,plotWithMask)

    MI = zeros(size(MeanAmp,1),1);
    MISurr = zeros(Nboot,size(MeanAmp,1));
    MXS = zeros(Nboot,1);
    
    if sum(sum(MeanAmp~=0))>0
        
        %% calc Modulation Index for original data
        MeanAmp = MeanAmp./repmat(sum(MeanAmp,2),1,size(MeanAmp,2));
        MeanAmp(isnan(MeanAmp)) = 0;
        MI = (log(nbBins)-(-sum(MeanAmp.*log(MeanAmp),2)))/log(nbBins);
        MI(isnan(MI)) = 0;

        if plotWithMask
            
            %% calc Modulation Index for surrogate data
            MeanAmpSurr = MeanAmpSurr./repmat(sum(MeanAmpSurr,3),1,1,size(MeanAmpSurr,3));
            MeanAmpSurr(isnan(MeanAmpSurr)) = 0;
            MISurr = (log(nbBins)-(-sum(MeanAmpSurr.*log(MeanAmpSurr),3)))/log(nbBins);
            MISurr(isnan(MISurr)) = 0;
            
            %% threshold MeanAmp values using Exteme Stat
            for fa = 1:size(MeanAmp,1)
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
            MXS = max(MISurr,[],2);
        end
    end
end