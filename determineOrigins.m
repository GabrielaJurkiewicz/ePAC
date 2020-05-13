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

function [CouplingOrigins] = determineOrigins(AverageSpectrum,SpectrumOfAverage,comodulogram,f,fA,fP,mask,w,Fs,epoch,fileID)
  
    CouplingOrigins = zeros(size(mask)); % 1 physiological, 0 epiphenomenal
    for fp = 1:size(mask,2)
        if sum(mask(:,fp)>0)>0
            ind = mask(:,fp)';
            ind = ind(2:end)-ind(1:end-1);
            start = find(ind==1);
            stop = find(ind==-1);
            if isempty(start)
                start = [1];
            else
                start = start+1;
            end
            if isempty(stop)
                stop = [size(mask,1)];
            else
                stop = stop;
            end
            if stop(1)<start(1)
                start = [1 start]; 
            end
            if start(end)>stop(end)
                stop = [stop size(mask,1)];
            end
            for ss=1:length(start)
                
                [~,idMaxPAC] = max(comodulogram(start(ss):stop(ss),fp));
                if isempty(idMaxPAC) || (idMaxPAC+start(ss)-1)==1
                    CouplingOrigins(start(ss):stop(ss),fp) = 0;
                    if (idMaxPAC+start(ss)-1)==1
                        txt = ['WARNING! Epoch ' num2str(epoch) ', fP = ' num2str(fP(fp)) ' Hz, fA = ' num2str(fA(1)) ' Hz: this coupling is too close to the lower amplitude frequency limit, thus it is labeled as Ambiguous. If you want to examine it, decrease the lower amplitude frequency limit.'];
                        fprintf(fileID,[txt '\n'])
                        disp(txt)
                    end
                else
                    wd = findWaveletFreqWidth(fA(idMaxPAC+start(ss)-1),w,Fs);
                    idfMAX = ((f>=(fA(idMaxPAC+start(ss)-1)-wd/2)) & (f<=(fA(idMaxPAC+start(ss)-1)+wd/2)));
                    idf = (f>=fA(start(ss)))&(f<=fA(stop(ss)));
                    ID = find(idfMAX|idf);
                    
                    AS = AverageSpectrum(ID,fp)/sum(AverageSpectrum(:,fp));
                    SA = SpectrumOfAverage(ID,fp)/sum(SpectrumOfAverage(:,fp));
                    tmpAS = AS(AS>SA);
                    SA = sum(SA(SA>AS));
                    AS = sum(tmpAS);

                    if (SA>=AS || SA>(0.999999*AS))
                        [~,idsa] = findpeaks(SpectrumOfAverage(ID,fp)');
                        [~,idM] = max(SpectrumOfAverage(ID,fp)');
                        if sum(idsa==idM)>0
                            idsa = idsa + ID(1) - 1;
                            [~,idmsa] = max(SpectrumOfAverage(idsa,fp));
                            idsa = idsa(idmsa);
                            if ~isempty(idsa)
                                if idfMAX(idsa)==1
                                    CouplingOrigins(start(ss):stop(ss),fp) = 1;
                                else
                                    CouplingOrigins(start(ss):stop(ss),fp) = 0;
                                end
                            else
                                CouplingOrigins(start(ss):stop(ss),fp) = 0;
                            end
                        else
                            CouplingOrigins(start(ss):stop(ss),fp) = 0;
                        end
                    else
                        [~,idas] = findpeaks(AverageSpectrum(ID,fp)');
                        [~,idM] = max(AverageSpectrum(ID,fp)');
                        if sum(idas==idM)>0
                            idas = idas + ID(1) - 1;
                            [~,idmas] = max(AverageSpectrum(idas,fp));
                            idas = idas(idmas);
                            if ~isempty(idas)
                                if idfMAX(idas)==1
                                    CouplingOrigins(start(ss):stop(ss),fp) = 1;
                                else
                                    CouplingOrigins(start(ss):stop(ss),fp) = 0;
                                end
                            else
                                CouplingOrigins(start(ss):stop(ss),fp) = 0;
                            end
                        else
                            CouplingOrigins(start(ss):stop(ss),fp) = 0;
                        end
                    end
                end
            end
        end
    end
    
    [B,L] = bwboundaries(mask,'noholes');
    for b = 1:length(B)
        centerFP = sum(sum(comodulogram.*(L==b).*repmat(fP,length(fA),1)))/sum(sum((L==b).*comodulogram));
        centerFA = sum(sum(comodulogram.*(L==b).*repmat(fA',1,length(fP))))/sum(sum((L==b).*comodulogram));
%         centerFP = round(mean(fP(sum(comodulogram==max(max(comodulogram.*(L==b))),1)>0)));
%         centerFA = round(mean(fA(sum(comodulogram==max(max(comodulogram.*(L==b))),2)>0)));
        originsB = sum(CouplingOrigins(L==b))>0;
        idBs = L(:,fP==round(centerFP));
        idBs = unique(idBs(idBs>0));
        idBs = setdiff(idBs,b);
        originsBs = [];
        centerFABs = [];
        if ~isempty(idBs)
            for ibs = 1:length(idBs)
                iifp = (fP==round(centerFP));
                iifa = L(:,fP==round(centerFP))==idBs(ibs);
%                 fAa = fA(iifa);
                originsBs = [originsBs sum(CouplingOrigins(iifa,iifp)>0)>0];
                centerFABs = [centerFABs sum(comodulogram(iifa,iifp)'.*fA(1,iifa))/sum(comodulogram(iifa,iifp))];
%                 centerFABs = [centerFABs fAa(comodulogram(iifa,iifp)==max(comodulogram(iifa,iifp)))];
            end
        end
        coef = 2;
        while round(centerFP*coef)<=fP(end)
            idB = L(:,fP==(round(centerFP*coef)));
            idB = unique(idB(idB>0));
            if ~isempty(idB)
                idB = setdiff(idB,b);
                if ~isempty(idB)
                    if ~isempty(originsBs)
                        for ib = 1:length(idB)
                            iifp = (fP==round(centerFP*coef));
                            iifa = L(:,fP==round(centerFP*coef))==idB(ib);
%                             fAa = fA(iifa);
                            centerfa = sum(comodulogram(iifa,iifp)'.*fA(1,iifa))/sum(comodulogram(iifa,iifp));
%                             centerfa = fAa(comodulogram(iifa,iifp)==max(comodulogram(iifa,iifp)));
                            if sum(abs(centerfa-centerFA)>abs(centerFABs-centerfa))==0
                                origBs = originsB;
                            else
                                [~,idfa] = min(abs(centerFABs-centerfa));
                                origBs = originsBs(idfa);
                            end
                            originsHarmonic = sum(CouplingOrigins(L==idB(ib)))>0;
                            fp1 = min(fP(sum(L==idB(ib),1)>0));
                            fp2 = max(fP(sum(L==idB(ib),1)>0));
                            fa1 = min(fA(sum(L==idB(ib),2)>0));
                            fa2 = max(fA(sum(L==idB(ib),2)>0));
                            if origBs==0 && originsHarmonic==1
                                CouplingOrigins(L==idB(ib)) = 0;
                                txt = ['WARNING! Epoch ' num2str(epoch) ', fP = <' num2str(fp1) ', ' num2str(fp2) '> Hz, fA = <' num2str(fa1) ', ' num2str(fa2) '> Hz: this coupling region appears to be a harmonic of an Ambiguous coupling, thus it is also labeled as Ambiguous.'];
                                fprintf(fileID,[txt '\n'])
                                disp(txt)
                            end
                        end
                    else
                        for ib = 1:length(idB)
                            originsHarmonic = sum(CouplingOrigins(L==idB(ib)))>0;
                            fp1 = min(fP(sum(L==idB(ib),1)>0));
                            fp2 = max(fP(sum(L==idB(ib),1)>0));
                            fa1 = min(fA(sum(L==idB(ib),2)>0));
                            fa2 = max(fA(sum(L==idB(ib),2)>0));
                            if originsB==0 && originsHarmonic==1
                                CouplingOrigins(L==idB(ib)) = 0;
                                txt = ['WARNING! Epoch ' num2str(epoch) ', fP = <' num2str(fp1) ', ' num2str(fp2) '> Hz, fA = <' num2str(fa1) ', ' num2str(fa2) '> Hz: this coupling region appears to be a harmonic of an Ambiguous coupling, thus it is also labeled as Ambiguous.'];
                                fprintf(fileID,[txt '\n'])
                                disp(txt)      
                            end
                        end
                    end
                end
            end
            coef = coef+1;
        end
    end
    
       
end