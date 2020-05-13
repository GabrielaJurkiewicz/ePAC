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
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%    Gabriela Jurkiewicz <gabriela.j.jurkiewicz@gmail.com>

function [Maps, MeanSignal, MeanFPsignal, MeanSpectrum, MIMaps, MIMessMaps, MIMeanFPsignal, num, numMI] = calcMeanMap(signal, MX, fPsignal, Fs, cut, cutMI, w, fAstart, fAend, fAstep, fA, fP_b, f, Nboot, plotWithMask, type)

    %% ---------------- PREPARE STRUCTURES ---------------------
    num   = 0;
    numMI = 0;
    Maps         = zeros(size(fA,2),floor(cut*Fs)*2+1);
    MeanSignal   = zeros(1,floor(cut*Fs)*2+1);
    MeanFPsignal = zeros(1,floor(cut*Fs)*2+1);
    MeanSpectrum = zeros(1,length(f));
    MIMaps       = zeros(size(fA,2),floor(cutMI*Fs)*2+1);
    MIMessMaps   = zeros(Nboot,size(fA,2),floor(cutMI*Fs)*2+1);
    MIMeanFPsignal = zeros(1,floor(cutMI*Fs)*2+1);
    
    %% ---------------- CALC MORLET WAVELET TF MAP ---------------------
    [map,~,~] = tf_cwt(hilbert(signal),fAstart,fAend+fAstep,Fs,w,fAstep,0); 
    
    %% ------------- ADDING MAPS ALIGNED TO EACH MAXIMA -------------- 
    mx = MX;
    while ~isempty(mx)
        
        idx  = mx(1,1);
        numMI = numMI + 1;
        mx(mx<(idx+floor(2*cutMI*Fs))) = [];
        
        MIMaps = MIMaps + map(:,idx-floor(cutMI*Fs):idx+floor(cutMI*Fs));
        MIMeanFPsignal  = MIMeanFPsignal + fPsignal(idx-floor(cutMI*Fs):idx+floor(cutMI*Fs));

        if plotWithMask
            for nb = 1:Nboot
                id = floor((rand()-1.0/2)*Fs/mean(fP_b));
                scale = (rand(1)-0.5)*0.1/0.5; % form -0.1 to 0.1
                scale = (1+scale);
                idxBoot = floor(idx*scale) + id;
                MAX = floor(size(map,2)*scale);
                while (idxBoot-floor(cutMI*Fs)<=0 || idxBoot+floor(cutMI*Fs) > MAX)
                    id = floor((rand()-1.0/2)*Fs/mean(fP_b));
                    scale = (rand(1)-0.5)*0.1/0.5; % form -0.1 to 0.1
                    scale = (1+scale);
                    idxBoot = floor(idx*scale) + id;
                    MAX = floor(size(map,2)*scale);
                end
                tmp = imresize(map,[size(map,1) floor(size(map,2)*scale)]);
                MIMessMaps(nb,:,:) = squeeze(MIMessMaps(nb,:,:)) + tmp(:,idxBoot-floor(cutMI*Fs):idxBoot+floor(cutMI*Fs));
            end
        end
    end
    
    %% ------------- ADDING MAPS ALIGNED TO nbCycles MAXIMA -------------- 
    mx = MX;
    while ~isempty(mx)
        
        idx  = mx(1,1);
        num = num + 1;
        mx(mx<(idx+floor(2*cut*Fs))) = [];
        
        s = signal(idx-floor(cut*Fs):idx+floor(cut*Fs));
        MeanSignal = MeanSignal + s;
        Maps = Maps + map(:,idx-floor(cut*Fs):idx+floor(cut*Fs));
        MeanFPsignal  = MeanFPsignal + fPsignal(idx-floor(cut*Fs):idx+floor(cut*Fs));
        [pxx,~] = periodogram(s,blackmanharris(length(s)),f,Fs);
        MeanSpectrum = MeanSpectrum + pxx;

    end

    %% ------------- MINIMAL NUMBER OF TRIALS CONDITION -------------------
    switch type
        case 'continuous'
            minNUM = 3;
        case 'event_related'
            minNUM = 1;
    end
    if (numMI < minNUM) || (num < 1)
        Maps           = zeros(size(fA,2),floor(cut*Fs)*2+1);
        MeanSignal     = zeros(1,floor(cut*Fs)*2+1);
        MeanFPsignal   = zeros(1,floor(cut*Fs)*2+1);
        MeanSpectrum   = zeros(1,length(f));
        MIMaps         = zeros(size(fA,2),floor(cutMI*Fs)*2+1);
        MIMessMaps     = zeros(Nboot,size(fA,2),floor(cutMI*Fs)*2+1);
        MIMeanFPsignal = zeros(1,floor(cutMI*Fs)*2+1);
        numMI          = 0;
        num            = 0;
    else
        if strcmp(type,'continuous')
            %% ------------- NORMALIZATION -------------------
            Maps         = Maps/num;
            MeanSignal   = MeanSignal/num;
            MeanFPsignal = MeanFPsignal/num;
            MeanSpectrum = MeanSpectrum/num;
            MIMaps       = MIMaps/numMI;
            if plotWithMask
                MIMessMaps = MIMessMaps/numMI;
            end
            MIMeanFPsignal = MIMeanFPsignal/numMI;
        end
    end

    
end
