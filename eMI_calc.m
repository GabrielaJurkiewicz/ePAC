%    Phase-amplitude coupling detection MATLAB plugin
%
%    Copyright (C) 2019 Gabriela Jurkiewicz
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%    Gabriela Jurkiewicz <gabriela.j.jurkiewicz@gmail.com>

function [] = eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,pPhaseCom,peMI,artifacts,plotWithMask,plotWithoutMask)

    %% -------------------------- RESTRICTIONS -------------------------
    if ((length(size(EEG.data))<3)&&(strcmp(type,'event_related')))
        disp('There are no epochs in this dataset. Compulsory change of "type" to continuous')
        type = 'continuous';
    end
    
    if (fAend>floor(EEG.srate/2))
        disp(['The highest amplitude frequency is the Nyquista frequency. Hence the compulsory change of "Frequency for amplitude stop" to ' num2str(floor(EEG.srate/2))])
        fAend = floor(EEG.srate/2);
    end
    
    if ~(plotWithMask || plotWithoutMask)
        disp('You checked not to plot any results. This choice is not optimal, thus the images with statistical mask will be saved.')
        plotWithMask = 1;
    end
    
    if ((fPstart-fPband/2)<1)
        disp(['The lowest frequency for phase with half of bandwidth subtracted is too close to 0 Hz. Compulsory change of "Frequency for phase start" to ' num2str(round(1+fPband/2))])
        fPstart = round(1+fPband/2);
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
    order = 4;
    nbCycles = 3;
    pSNR = 95;
    margines = w/fAstart;
    filtDistortion = order*2/EEG.srate;
    if filtDistortion > margines
        margines = filtDistortion;
    end
    fP = fPstart:fPstep:fPend;
    fP_bins = zeros(2,length(fP));
    fP_bins(1,:) = fP-fPband/2;
    fP_bins(2,:) = fP+fPband/2;
    [status,message,~] = mkdir(dirOut);
    [ID] = findPreviousVersion(dirOut,'eMI');
    
    %% -------------------- SAVE CONFIGURATION FILE -------------------------
    disp('----------------------- eMI: saving configuration file ------------------------------------------------------')
    config.EEGname     = EEG.filename;
    config.EEGpath     = EEG.filepath;
    config.EEGdataHash = DataHash(EEG.data,struct('Input','array'));
    config.Fs          = EEG.srate;
    config.chan_fP     = chan_fP;
    config.chan_fA     = chan_fA;
    config.epochs      = epochs;
    config.fPstart     = fPstart;
    config.fPend       = fPend;
    config.fPstep      = fPstep;
    config.fPband      = fPband;
    config.fAstart     = fAstart;
    config.fAend       = fAend;
    config.fAstep      = fAstep;
    config.dirOut      = dirOut;
    config.type        = type;
    config.nbCycles    = nbCycles;
    config.nbBins      = nbBins;
    config.w           = w;
    config.Nboot       = Nboot;
    config.pPhaseCom   = pPhaseCom;
    config.peMI        = peMI;
    config.artifacts   = artifacts;
    config.plotWithMask = plotWithMask;
    config.plotWithoutMask = plotWithoutMask;
    save([dirOut 'v' ID '_eMI_config.mat'],'config')

    %% ------------ FIND LOW-FREQ OSCILLATION AND IT'S MAXIMA ---------------
    disp('----------------------- eMI: preparing low-frequency oscillation --------------------------------------------')
    [Maxes, LowFreq] = findLowFreqMaxes(EEG,chan_fP,fP_bins,Nboot,order,pSNR);
    HighFreqSignal   = EEG.data(chan_fA,:,:);
    
    %% ---------------------------- MAIN PART -------------------------------
    disp('----------------------- eMI: calculations -------------------------------------------------------------------')
    if strcmp(type,'continuous')
        run_continuous_eMI(EEG,Epochs,HighFreqSignal,Maxes,LowFreq,fP_bins,fP,fAstart,fAend,fAstep,margines,dirOut,nbCycles,w,nbBins,Nboot,pPhaseCom,peMI,ID,plotWithMask,plotWithoutMask)
    elseif strcmp(type,'event_related')
        run_event_related_eMI(EEG,Epochs,HighFreqSignal,Maxes,LowFreq,fP_bins,fP,fAstart,fAend,fAstep,margines,dirOut,nbCycles,w,nbBins,Nboot,pPhaseCom,peMI,ID,plotWithMask,plotWithoutMask)
    end    
    
end
