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

function [] = run_event_related_dPAC(EEG,SignalA,SignalP,Epochs,fP_bins,fP,fA_bins,fA,dirOut,Nboot,pdPAC,ID,plotWithMask,plotWithoutMask)
    
    
    %% -------------------- PREPARE VARIABLES ----------------------------- 
    Fs           = EEG.srate;
    comodulogram = zeros(length(fP), length(fA));
    surrogates   = zeros(Nboot, length(fP), length(fA));
    pval         = zeros(length(fP), length(fA)); 

    results = matfile([dirOut 'v' ID '_dPAC_results.mat'],'Writable',true);
    results.fA = fA;
    results.fP = fP;
    results.fP_bins = fP_bins;
    results.fA_bins = fA_bins;       
        
    progressbar('on')
    
    %% -------------------- EXTRACT PHASE AND AMPLITUDE SIGNAL -----------------------------
    Amplitude = zeros(length(fA),length(Epochs),size(SignalA,2));
    Phase = zeros(length(fP),length(Epochs),size(SignalP,2));
    
    for epoch = 1:length(Epochs)
        for i = 1:length(fP)
            low = eegfilt(squeeze(SignalP(1,:,Epochs(epoch))),Fs,fP_bins(1,i),fP_bins(2,i)); 
            Phase(i,epoch,:) = angle(hilbert(low));
        end
        for j = 1:length(fA)
            high = eegfilt(squeeze(SignalA(1,:,Epochs(epoch))),Fs,fA_bins(1,j),fA_bins(2,j));  
            Amplitude(j,epoch,:) = abs(hilbert(high));                 
        end
        progressbar(epoch/length(Epochs)*100/10)
    end
    
    %% -------------------- CALCULATE MODULATION INDEX -----------------------------
    %   adapted from the code written by Tolga Ozkurt (Ozkurt and Schnitzler, 2011)
    for i = 1:length(fP)
        tmpP = Phase(i,:,Fs:size(Phase,3)-Fs-1);
        for j = 1:length(fA)
            tmpA = Amplitude(j,:,Fs:size(Phase,3)-Fs-1);
            comodulogram(i,j) = abs(mean(tmpA(:).*exp(1i*tmpP(:)))) / sqrt(mean(tmpA(:).*tmpA(:)));
            if plotWithMask
                for nb = 1:Nboot
                    PhaseSurr = zeros(length(Epochs)*size(tmpP,3),1);
                    for epoch = 1:length(Epochs)                  
                        Surrogate = makeSurrogateNoiseFilt(Phase(i,epoch,:),fP_bins(1,i),fP_bins(2,i),Fs,1);
                        PhaseSurr((epoch-1)*size(tmpP,3)+1:epoch*size(tmpP,3),1) = angle(hilbert(Surrogate(1,Fs:size(Phase,3)-Fs-1)));
                    end
                    surrogates(nb,i,j) = abs(mean(tmpA(:).*exp(1i*PhaseSurr))) / sqrt(mean(tmpA(:).*tmpA(:)));
                end 
            end
            progressbar(10+9*((i-1)*length(fA)+j)/(length(fP)*length(fA))*100/10)
        end
    end
    
    %% -------------------- STATISTICS -----------------------------
    if plotWithMask
        MX = squeeze(max(max(surrogates,[],2),[],3));
        thresh = prctile(MX,pdPAC);
        comodulogramStat = comodulogram;
        comodulogramStat(comodulogramStat<thresh) = 0;
        for i = 1:length(fP)
            for j = 1:length(fA)
                pval(i,j) = sum(MX>=comodulogram(i,j))/Nboot;
            end
        end
        comodulogramStat = comodulogramStat';
        pval             = pval';
    end
    comodulogram = comodulogram';

    %% -------------------- SAVE AND PLOT RESULTS + CLEAN -----------------------------
    results.comodulogram = comodulogram;
    if plotWithMask
        results.comodulogramStat = comodulogramStat;
        results.pval = pval;
    end
    matfile([dirOut 'v' ID '_dPAC_results.mat'],'Writable',false);

    clear pval results

    if plotWithoutMask
        plotFigure(fP,fA,comodulogram,[],'dPAC',dirOut,['v' ID '_dPAC_comodulogram'])
    end
    if plotWithMask
        plotFigure(fP,fA,comodulogramStat,[],'dPAC',dirOut,['v' ID '_dPAC_comodulogramStat'])
    end
    progressbar('off')
    
end