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

function [] = run_continuous_MI(EEG,SignalA,SignalP,Epochs,fP_bins,fP,fA_bins,fA,dirOut,Nboot,pMI,ID,phaseBins,plotWithMask,plotWithoutMask)
    
    Fs = EEG.srate;
    for epoch = 1:length(Epochs)
        
        disp(['Epoch ' num2str(Epochs(epoch)) ' / ' num2str(EEG.trials)])
        progressbar('on')
        
        %% -------------------- PREPARE VARIABLES -----------------------------
        comodulogram = zeros(length(fP), length(fA));
        surrogates   = zeros(Nboot, length(fP), length(fA));
        pval         = zeros(length(fP), length(fA)); 
        signalA      = squeeze(SignalA(1,:,Epochs(epoch)));
        signalP      = squeeze(SignalP(1,:,Epochs(epoch)));
        
        results = matfile([dirOut 'v' ID '_MI_epoch' num2str(Epochs(epoch)) '_results.mat'],'Writable',true);
        results.fA = fA;
        results.fP = fP;
        results.fP_bins = fP_bins;
        results.fA_bins = fA_bins;       
        
        %% -------------------- CALCULATE MODULATION INDEX -----------------------------
        % based on Tort et al., J Neurophysiol 2010
        for i = 1:length(fP)

            theta = eegfilt(signalP,Fs,fP_bins(1,i),fP_bins(2,i));
            phase = angle(hilbert(theta));                              

            for j=1:length(fA)

                gamma = eegfilt(signalA,Fs,fA_bins(1,j),fA_bins(2,j));   
                amp   = abs(hilbert(gamma));                             
                comodulogram(i,j) = MI(phase, amp, phaseBins);
                
                if plotWithMask
                    Surrogate = makeSurrogateNoiseFilt(signalP,fP_bins(1,i),fP_bins(2,i),Fs,Nboot);
                    for nb = 1:Nboot
                        phaseS = angle(hilbert(Surrogate(nb,:)));
                        surrogates(nb,i,j) = MI(phaseS, amp, phaseBins);
                    end
                end
                progressbar(( (i-1)/length(fP) + j/(length(fA)*length(fP)) ) * 100)
                
            end
            progressbar(i/length(fP)*100)
            
        end
        
        %% -------------------- STATISTICS -----------------------------
        if plotWithMask
            MX = squeeze(max(max(surrogates,[],2),[],3));
            thresh = prctile(MX,pMI);
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
        matfile([dirOut 'v' ID '_MI_epoch' num2str(Epochs(epoch)) '_results.mat'],'Writable',false);
        
        clear pval results
        
        if plotWithoutMask
            plotFigure(fP,fA,comodulogram,[],'MI',dirOut,['v' ID '_MI_epoch' num2str(Epochs(epoch)) '_comodulogram'])
        end
        if plotWithMask
            plotFigure(fP,fA,comodulogramStat,[],'MI',dirOut,['v' ID '_MI_epoch' num2str(Epochs(epoch)) '_comodulogramStat'])
        end
        progressbar('off')
        
    end

    if length(Epochs)>1
        plotContinuousSummary(EEG,fP,fA,dirOut,ID,(100-pmMVL)/100,'MI',plotWithMask,plotWithoutMask)
    end

end
