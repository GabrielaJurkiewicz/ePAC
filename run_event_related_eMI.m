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

function [] = run_event_related_eMI(EEG,Epochs,HighFreqSignal,Maxes,LowFreq,fP_bins,fP,fAstart,fAend,fAstep,margines,dirOut,nbCycles,w,nbBins,Nboot,pPhaseCom,peMI,ID,plotWithMask,plotWithoutMask)

    progressbar('on')
    fA = fAstart:fAstep:fAend;
    fSpectrum = fAstart:1:fAend;
    phaseBins = zeros(1,nbBins);
    dPhi      = 2*pi/nbBins;
    for b = 1:nbBins 
        phaseBins(b) = -pi+(b-1)*dPhi; 
    end
     
    comodulogram = zeros(length(fA),size(fP_bins,2));
    comodulogramSurr = zeros(Nboot,length(fA),size(fP_bins,2));
    phaseComodulogram = zeros(length(fA),nbBins,size(fP_bins,2));
    averageSpectrum = zeros(length(fSpectrum),length(fP));
    spectrumOfAverage = zeros(length(fSpectrum),length(fP));
        
    results = matfile([dirOut 'v' ID '_eMI_results.mat'],'Writable',true);
    results.fA = fA;
    results.fP = fP;
    results.fP_bins = fP_bins;
    results.fSpectrum = fSpectrum;
    
    for i = 1:size(fP_bins,2)

        %% ---------------- INITIATE VARIABLES  ---------------------  
        fP_b  = fP_bins(:,i)';
        cut   = nbCycles/(2*fP(1,i));
        cutMI = 1/(2*fP(1,i));
        Time  = linspace(-cut,cut,floor(cut*EEG.srate)*2+1);
        Map          = zeros(size(fA,2),floor(cut*EEG.srate)*2+1);
        MeanSignal   = zeros(1,floor(cut*EEG.srate)*2+1);
        MeanLowFreq  = zeros(1,floor(cut*EEG.srate)*2+1);
        MeanSpectrum = zeros(1,length(fSpectrum));
        MeanAmp      = zeros(length(fA),nbBins);
        MeanAmpSurr  = zeros(Nboot,length(fA),nbBins);
        
        Num = 0;
        NumMI = 0;
        for epoch = 1:length(Epochs)
            
            %% ---------------- ADJUST DATA TO SPECIFIC EPOCH AND FP ---------------------
            signal   = squeeze(HighFreqSignal(1,:,Epochs(epoch)));
            fPsignal = LowFreq{i}{1,Epochs(epoch)};
            mxs      = Maxes{i}{1,Epochs(epoch)};
            mx       = mxs((mxs>(floor((cut+margines)*EEG.srate)))&(mxs<(floor(size(signal,2))-floor((cut+margines)*EEG.srate))));
            
            %% ---------------- CALCULATE TF MAPS & ASSIGN PHASE BINS -------------------
            [map, meanSignal, meanLowFreq, meanSpectrum, MImap, MImessMap, MImeanLowFreq, num, numMI] = calcMeanMap(signal, mx, fPsignal, EEG.srate, cut, cutMI, w, fAstart, fAend, fAstep, fA, fP_b, fSpectrum, Nboot, plotWithMask, 'event_related');
            Num = Num + num;
            NumMI = NumMI + numMI;
            Map = Map + map;
            MeanSignal = MeanSignal + meanSignal;
            MeanLowFreq = MeanLowFreq + meanLowFreq;
            MeanSpectrum = MeanSpectrum + meanSpectrum;
            clear map meanSignal meanLowFreq meanSpectrum
            
            [meanAmp,meanAmpSurr] = assignPhaseBins(MImap,MImessMap,MImeanLowFreq,nbBins,phaseBins,Nboot,plotWithMask);
            MeanAmp = MeanAmp + meanAmp;
            MeanAmpSurr = MeanAmpSurr + meanAmpSurr;
            clear MImap MImessMap MImeanLowFreq meanAmp meanAmpSurr 
            
            progressbar(((i-1)*length(Epochs)+epoch)/(size(fP_bins,2)*length(Epochs))*100)
            
        end 
        %% -------------------- NORMALIZE -----------------------------
        if NumMI < 3 
            Map          = zeros(size(fA,2),floor(cut*EEG.srate)*2+1);
            MeanSignal   = zeros(1,floor(cut*EEG.srate)*2+1);
            MeanLowFreq  = zeros(1,floor(cut*EEG.srate)*2+1);
            MeanSpectrum = zeros(1,length(fSpectrum));
            MeanAmp      = zeros(length(fA),nbBins);
            MeanAmpSurr  = zeros(Nboot,length(fA),nbBins);
        else
            Map          = Map/Num;
            MeanSignal   = MeanSignal/Num;
            MeanLowFreq  = MeanLowFreq/Num;
            MeanSpectrum = MeanSpectrum/Num;
            MeanAmp      = MeanAmp/NumMI;
            MeanAmpSurr  = MeanAmpSurr/NumMI;
        end
        
        %% -------------------- CALCULATE SPECTRUM OF MEAN SIGNAL -----------------------------
        averageSpectrum(:,i) = MeanSpectrum;
        [pxx,~] = periodogram(MeanSignal,blackmanharris(length(MeanSignal)),fSpectrum,EEG.srate);
        spectrumOfAverage(:,i) = pxx;
        
        %% -------------------- CALCULATE MI -----------------------------
        [MI,MISurr,PC] = calc_ModulationIndex_ER(MeanAmp,MeanAmpSurr,nbBins,Nboot,pPhaseCom,plotWithMask);
        comodulogramSurr(:,i) = MISurr;
        comodulogram(:,i) = MI;
        phaseComodulogram(:,:,i) = PC;  
        clear MI MISurr PC
        
        %% -------------------- SAVE AND PLOT RESULTS + CLEAN -----------------------------
        eval(['results' num2str(fP(1,i)) 'Hz.Time = Time;'])
        eval(['results' num2str(fP(1,i)) 'Hz.Map = Map;'])
        eval(['results' num2str(fP(1,i)) 'Hz.MeanSignal = MeanSignal;'])
        eval(['results' num2str(fP(1,i)) 'Hz.MeanLowFreq = MeanLowFreq;'])
        eval(['results.results' num2str(fP(1,i)) 'Hz = results' num2str(fP(1,i)) 'Hz;'])
        eval(['clear results' num2str(fP(1,i)) 'Hz'])
        clear Map MeanSignal MeanLowFreq

        progressbar(i/size(fP_bins,2)*100)
    end
    
    progressbar('off')
    disp('----------------------- eMI: ploting and saving results -----------------------------------------------------')
        
    %% ---------------- eMI STATISTIC --------------------------------------
    if plotWithMask
        MX  = max(max(comodulogramSurr,[],2),[],3);
        clear comodulogramSurr
        tresh = prctile(MX,peMI);
        mask = (comodulogram>=tresh);
        mask = mask & squeeze(sum(phaseComodulogram>0,2)>0);

        comodulogramStat = zeros(size(comodulogram));
        comodulogramStat(mask) = comodulogram(mask);
        
        pval = zeros(size(comodulogram));
        phaseComodulogramStat = zeros(size(phaseComodulogram));
        for fp = 1:length(fP)
            for fa = 1:length(fA)
                pval(fa,fp) = sum(MX>=comodulogram(fa,fp))/Nboot;
                if mask(fa,fp)==1
                    phaseComodulogramStat(fa,:,fp) = phaseComodulogram(fa,:,fp);
                end
            end
        end
        
        fileID = fopen([dirOut 'v' ID '_eMI_WARNINGS.txt'],'w');
        couplingOrigins = determineOrigins(averageSpectrum,spectrumOfAverage,comodulogram,fSpectrum,fA,fP,mask,w,EEG.srate,Epochs(epoch),fileID);
        fclose(fileID);
    end
        
    %% ---------------- SAVE AND PLOT RESULTS + CLEAN --------------------------------------
    if plotWithoutMask
        plotFigure(fP,fA,comodulogram,[],'eMI',dirOut,['v' ID '_eMI_comodulogram'])
%         plotFigurePhase(fP,fA,phaseBins,phaseComodulogram,[],'eMI',dirOut,['v' ID '_eMI_phaseComodulogram']) 
    end
    results.phaseComodulogram = phaseComodulogram;
    clear phaseComodulogram
    results.averageSpectrum = averageSpectrum;
    clear averageSpectrum
    results.spectrumOfAverage = spectrumOfAverage;
    clear spectrumOfAverage
    results.comodulogram = comodulogram;
    clear comodulogram
        
    if plotWithMask
        plotFigure(fP,fA,comodulogramStat,couplingOrigins,'eMI',dirOut,['v' ID '_eMI_comodulogramStat'])
%         plotFigurePhase(fP,fA,phaseBins,phaseComodulogramStat,couplingOrigins,'eMI',dirOut,['v' ID '_eMI_phaseComodulogramStat']) 
        plotPolarHist(fP,fA,phaseBins,comodulogramStat,phaseComodulogramStat,couplingOrigins,nbBins,dirOut,['v' ID '_eMI_phaseHistogramStat']) 
            
        results.phaseComodulogramStat = phaseComodulogramStat;
        clear phaseComodulogramStat
        results.comodulogramStat = comodulogramStat;
        clear comodulogramStat
        results.pval = pval;
        clear pval
        results.couplingOrigins = couplingOrigins;
        clear couplingOrigins
    end
    
    results.phaseBins = phaseBins;
    clear phaseBins
    matfile([dirOut 'v' ID '_eMI_results.mat'],'Writable',false);
    clear results

    load([dirOut 'v' ID '_eMI_results.mat']);
    for i = 1:size(fP_bins,2)
        cut   = nbCycles/(2*fP(1,i));
        Time  = linspace(-cut,cut,floor(cut*EEG.srate)*2+1);
        eval(['Time = results' num2str(fP(1,i)) 'Hz.Time;'])
        eval(['Map = results' num2str(fP(1,i)) 'Hz.Map;'])
        eval(['MeanSignal = results' num2str(fP(1,i)) 'Hz.MeanSignal;'])
        eval(['MeanLowFreq = results' num2str(fP(1,i)) 'Hz.MeanLowFreq;'])
        if plotWithMask
            comod = comodulogramStat(:,i);
            SA = spectrumOfAverage(:,i);
            AS = averageSpectrum(:,i);
            CO = couplingOrigins(:,i);
            PC = squeeze(phaseComodulogramStat(:,:,i));
        else
            comod = []; SA = []; AS = []; CO = []; PC = zeros(length(fA),nbBins);
        end
        plot_maps_together(Time,fA,fSpectrum,Map,MeanSignal,MeanLowFreq,CO,PC,[phaseBins pi],SA,AS,comod,w,EEG.srate,...
                               [dirOut 'v' ID '_eMI_results' num2str(fP(1,i)) 'Hz'],cut,plotWithMask)
    end
    clear Time Map MeanSignal MeanLowFreq comodulogram comodulogramStat phaseComodulogramStat phaseComodulogram phaseBins spectrumOfAverage averageSpectrum couplingOrigins

end