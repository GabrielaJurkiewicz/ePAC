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

function [Maxes, LowFreq] = findLowFreqMaxes(EEG,chan_fP,fP_bins,Nboot,order,pSNR)

    progressbar('on')

    % create reference 'signal to noise ratio' of spectrum
    for nb = 1:Nboot
        sigNoise = pinknoise(1,size(EEG.data,2));
        if size(EEG.data,2) < 2*EEG.srate
            N = 2*EEG.srate-size(EEG.data,2);
            sigNoise = [zeros(1,ceil(N/2)) sigNoise zeros(1,ceil(N/2))];
        end
        [pxxN,f] = pwelch(sigNoise,hamming(2*EEG.srate),EEG.srate,(1:0.5:EEG.srate/2),EEG.srate);
        [pxxInterp] = fitBackgroundSpectrum(pxxN,f);
        SNRnoise  = pxxN./pxxInterp;
        if nb == 1
           SNRDistribution = zeros(Nboot,length(f)); 
        end
        SNRDistribution(nb,:) = SNRnoise;
        progressbar(nb/Nboot*30)
    end
    fPidx = (f>=fP_bins(1,1))&(f<fP_bins(2,end));
    snrd = SNRDistribution(:,fPidx);
    SNRthresh = prctile(snrd(:),pSNR);

    % check which low-frequency is above SNR threshold - is significant
    significantLowFreq = zeros(size(fP_bins,2),EEG.trials);
    for trial = 1:EEG.trials
        sig = double(EEG.data(chan_fP,:,trial));
        if size(EEG.data,2) < 2*EEG.srate
            N = 2*EEG.srate-size(EEG.data,2);
            sig = [zeros(1,ceil(N/2)) sig zeros(1,ceil(N/2))];
        end
        [pxx,f] = pwelch(sig,hamming(2*EEG.srate),EEG.srate,(1:0.5:EEG.srate/2),EEG.srate);
        [pxxInterp] = fitBackgroundSpectrum(pxx,f);
        SNR = pxx./pxxInterp;
        for Idx = 1:size(fP_bins,2)
            fP_b = fP_bins(:,Idx)';
            fPidx = (f>=fP_b(1))&(f<fP_b(2));
            if sum(SNR(fPidx)>SNRthresh)>0
                significantLowFreq(Idx,trial) = 1;
            end
            progressbar(30+((trial-1)*size(fP_bins,2)+Idx)/(size(fP_bins,2)*EEG.trials)*30)
        end
        progressbar(30+trial/EEG.trials*30)
    end
    
    % find significant-low-freq oscillations and their maxima
    Maxes = {};
    LowFreq = {};
    for Idx = 1:size(fP_bins,2)
        
        fP_b = fP_bins(:,Idx)';
        Fs   = EEG.srate; % Sampling Frequency
        N    = order;     % Order
        Fc1  = fP_b(1);   % First Cutoff Frequency
        Fc2  = fP_b(2);   % Second Cutoff Frequency
        h    = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
        Hd   = design(h, 'butter');
        Mx = {};
        LF = {};
        
        for trial = 1:EEG.trials
            if significantLowFreq(Idx,trial)==1
                sig = double(EEG.data(chan_fP,:,trial));
                sigF = filtfilt(Hd.sosMatrix,Hd.ScaleValues,sig);
                sigF = sigF-mean(sigF);
                [~,~,~,prom] = findpeaks(double(sigF));
                [~,IdxMax] = findpeaks(double(sigF),'MINPEAKPROMINENCE', median(prom)*0.05);
                if isempty(IdxMax)
                    IdxMax = [];
                end
            else
                sigF   = [];
                IdxMax = [];
            end
            Mx{trial} = IdxMax;
            LF{trial} = sigF;
            progressbar(60+((Idx-1)*EEG.trials+trial)/(size(fP_bins,2)*EEG.trials)*40)
        end
        
        Maxes{Idx}   = Mx;
        LowFreq{Idx} = LF;
        progressbar(60+Idx/size(fP_bins,2)*40)
        
    end
    progressbar('off')

end