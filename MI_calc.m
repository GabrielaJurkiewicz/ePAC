function [] = MI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,nbBins,Nboot,pMI,artifacts,plotWithMask,plotWithoutMask)

    %% -------------------------- RESTRICTIONS -------------------------
    if ((length(size(EEG.data))<3)&&(strcmp(type,'event_related')))
        disp('There are no epochs in this dataset. Compulsory change of "type" to continuous')
        type = 'continuous';
    end
    
    if (fAfiltband < 2*fPend)
        disp(['The amplitude-frequency filter width must be at least 2 times bigger than the highest phase-frequency. It is so, because the spectrum of coupled signal containes side peaks (fA-fP,fA+fP) around coupled high-frequency. Hence the compulsory change of "Frequency for amplitude filtration band width:" to ' num2str(2*fPend)])
        fAfiltband = 2*fPend;
    end
    
    if (fAstart < fPend)
        disp(['The minimal amplitude-frequency is equal to maximal phase-frequency. It is so, because the bandwidth of phase-frequency filter is constant for all amplitude-frequencies and it must be 2 times bigger than the highest phase-frequency. Hence the compulsory change of "Frequency for amplitude start:" to ' num2str(fPend)])
        fAstart = fPend;
    end

    if (fAend>floor(EEG.srate/2/1.15-fAfiltband/2))
        disp(['The highest amplitude frequency is the Nyquista frequency divided by 1.15 (accounts for transision zone), with a half of filtration bandwidth subtracted. Hence the compulsory change of "Frequency for amplitude stop" to ' num2str(floor(EEG.srate/2/1.15-fAfiltband/2))])
        fAend = floor(EEG.srate/2/1.15-fAfiltband/2);
    end
    
    if ~(plotWithMask || plotWithoutMask)
        disp(['You checked not to plot any results. This choice is not optimal, thus the images with statistical mask will be saved.'])
        plotWithMask = 1;
    end
    
    if ((fPstart-fPfiltband/2)<1)
        disp(['The lowest frequency for phase with half of bandwidth subtracted is too close to 0 Hz. Compulsory change of "Frequency for phase start" to ' num2str(round(1+fPfiltband/2))])
        fPstart = round(1+fPfiltband/2);
    end
    
    %% -------------------------- SELECT EPOCHS -------------------------
    idx = zeros(1,EEG.trials);
    if artifacts
        if ~isempty(EEG.reject.rejjpE)
            idx = idx | (sum(EEG.reject.rejjpE([chan_fP,chan_fA],:),1)>0);
        end
        if ~isempty(EEG.reject.rejkurtE)
            idx = idx | (sum(EEG.reject.rejkurtE([chan_fP,chan_fA],:),1)>0);
        end
        if ~isempty(EEG.reject.rejmanualE)
            idx = idx | (sum(EEG.reject.rejmanualE([chan_fP,chan_fA],:),1)>0);
        end
        if ~isempty(EEG.reject.rejthreshE)
            idx = idx | (sum(EEG.reject.rejthreshE([chan_fP,chan_fA],:),1)>0);
        end
        if ~isempty(EEG.reject.rejconstE)
            idx = idx | (sum(EEG.reject.rejconstE([chan_fP,chan_fA],:),1)>0);
        end
        if ~isempty(EEG.reject.rejfreqE)
            idx = idx | (sum(EEG.reject.rejfreqE([chan_fP,chan_fA],:),1)>0);
        end
    end
    
    idxUser = ones(1,EEG.trials);
    idxUser(epochs) = 0;
    idx = idx | idxUser;
    Epochs = find(~(idx));
                        
    %% -------------------------- PREPARE VARIABLES -------------------------
    fP = fPstart:fPstep:fPend;
    fP_bins = zeros(2,length(fP));
    fP_bins(1,:) = fP-fPfiltband/2;
    fP_bins(2,:) = fP+fPfiltband/2;

    fA = fAstart:fAstep:fAend;
    fA_bins = zeros(2,length(fA));
    fA_bins(1,:) = fA-fAfiltband/2;
    fA_bins(2,:) = fA+fAfiltband/2;
    
    % this variable will get the beginning (not the center) of each phase bin (in rads)
    phaseBins = zeros(1,nbBins); 
    for j = 1:nbBins 
        phaseBins(j) = -pi+(j-1)*2*pi/nbBins; 
    end
    
    signalA = double(EEG.data(chan_fA,:,:));
    signalP = double(EEG.data(chan_fP,:,:));
    
    [status,message,~] = mkdir(dirOut);
    [ID] = findPreviousVersion(dirOut,'MI');
    
    %% -------------------- SAVE CONFIGURATION FILE -------------------------
    config.EEGname     = EEG.filename;
    config.EEGpath     = EEG.filepath;
    config.EEGdataHash = DataHash(EEG.data,struct('Input','array'));
    config.chan_fP     = chan_fP;
    config.chan_fA     = chan_fA;
    config.epochs      = epochs;
    config.fPstart     = fPstart;
    config.fPend       = fPend;
    config.fPstep      = fPstep;
    config.fPfiltband  = fPfiltband;
    config.fAstart     = fAstart;
    config.fAend       = fAend;
    config.fAstep      = fAstep;
    config.fAfiltband  = fAfiltband;
    config.dirOut      = dirOut;
    config.type        = type;
    config.nbBins      = nbBins;
    config.Nboot       = Nboot;
    config.pMI         = pMI;
    config.artifacts   = artifacts;
    config.plotWithMask = plotWithMask;
    config.plotWithoutMask = plotWithoutMask;
    save([dirOut 'v' ID '_MI_config.mat'],'config')
    
    %% ---------------------------- MAIN PART -------------------------------
    disp('-----------------------Modulation Index algorithm calculations----------------------')
    if strcmp(type,'continuous')
        run_continuous_MI(EEG,signalA,signalP,Epochs,fP_bins,fP,fA_bins,fA,dirOut,Nboot,pMI,ID,phaseBins,plotWithMask,plotWithoutMask)
    elseif strcmp(type,'event_related')
        run_event_related_MI(EEG,signalA,signalP,Epochs,fP_bins,fP,fA_bins,fA,dirOut,Nboot,pMI,ID,phaseBins,plotWithMask,plotWithoutMask)
    end
    

end
