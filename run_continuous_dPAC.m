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

function [] = run_continuous_dPAC(EEG,SignalA,SignalP,Epochs,fP_bins,fP,fA_bins,fA,dirOut,Nboot,pdPAC,ID,plotWithMask,plotWithoutMask)
    
    Fs = EEG.srate;
    for epoch = 1:length(Epochs)
        
        disp(['Epoch ' num2str(Epochs(epoch)) ' / ' num2str(EEG.trials)])
        progressbar('on')
        
        %% -------------------- PREPARE VARIABLES -----------------------------
        comodulogram = zeros(length(fP), length(fA));
        surrogates   = zeros(Nboot, length(fP), length(fA));
        pval         = zeros(length(fP), length(fA)); 
        signalA      = squeeze(SignalA(1,:,epoch));
        signalP      = squeeze(SignalP(1,:,epoch));
        
        results = matfile([dirOut 'v' ID '_dPAC_epoch' num2str(Epochs(epoch)) '_results.mat'],'Writable',true);
        results.fA = fA;
        results.fP = fP;
        results.fP_bins = fP_bins;
        results.fA_bins = fA_bins;       
        
        %% -------------------- CALCULATE DIRECT PAC ESTIMATE -----------------------------
        %   adapted from the code written by Tolga Ozkurt (Ozkurt and Schnitzler, 2011)
        for i = 1:length(fP)

            theta = eegfilt(signalP,Fs,fP_bins(1,i),fP_bins(2,i));       %Compute the low freq signal.
            theta = theta(Fs:length(signalP)-Fs-1);                      %Drop the first and last second.
            phase = angle(hilbert(theta));                               %Compute the low freq phase.

            for j=1:length(fA)

                gamma = eegfilt(signalA,Fs,fA_bins(1,j),fA_bins(2,j));   %Compute the high freq signal.
                gamma = gamma(Fs:length(signalA)-Fs-1);                  %Drop the first and last second.
                amp   = abs(hilbert(gamma));                             %Compute the high freq amplitude.
                comodulogram(i,j) = abs(mean(amp.*exp(1i*phase))) / sqrt(mean(amp.*amp));
                
                if plotWithMask
                    Surrogate = makeSurrogateNoiseFilt(signalP,fP_bins(1,i),fP_bins(2,i),Fs,Nboot);
                    for nb = 1:Nboot
                        tmp = Surrogate(nb,:);
                        phaseS = angle(hilbert(tmp(Fs:length(signalP)-Fs-1)));
                        surrogates(nb,i,j) = abs(mean(amp.*exp(1i*phaseS))) / sqrt(mean(amp.*amp));
                    end
                end
                progressbar(( (i-1)/length(fP) + j/(length(fA)*length(fP)) ) * 100)
                
            end
            progressbar(i/length(fP)*100)
            
        end

        %% -------------------- STATISTICS -----------------------------
        if plotWithMask
            MX = squeeze(max(max(surrogates,[],2),[],3));
            thresh = prctile(MX,pdPAC);
            comodulogramStat = comodulogram;
            comodulogramStat(comodulogramStat<thresh) = 0;
            for i=1:length(fP)
                for j=1:length(fA)
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
            results.pval = pval;
            results.comodulogramStat = comodulogramStat;
        end
        matfile([dirOut 'v' ID '_dPAC_epoch' num2str(Epochs(epoch)) '_results.mat'],'Writable',false);
        
        clear pval results
        
        if plotWithoutMask
            plotFigure(fP,fA,comodulogram,[],'dPAC',dirOut,['v' ID '_dPAC_epoch' num2str(Epochs(epoch)) '_comodulogram'])
        end
        if plotWithMask
            plotFigure(fP,fA,comodulogramStat,[],'dPAC',dirOut,['v' ID '_dPAC_epoch' num2str(Epochs(epoch)) '_comodulogramStat'])
        end
        progressbar('off')
        
    end

    if length(Epochs)>1
        plotContinuousSummary(EEG,fP,fA,dirOut,ID,(100-pmMVL)/100,'dPAC',plotWithMask,plotWithoutMask)
    end

end
