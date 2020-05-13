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

function [] = plotContinuousSummary(EEG,fP,fA,dirOut,ID,alpha,methodName,plotWithMask,plotWithoutMask)

    fontSize = 12;
    fontlist = listfonts();
    if sum(strcmp(fontlist,'Arial'))>0
        fontNm = 'Arial';
    else
        if sum(strcmp(fontlist,'Liberation Sans'))>0
            fontNm = 'Liberation Sans';
        else
            if sum(strcmp(fontlist,'Helvetica'))>0
                fontNm = 'Helvetica';
            end
        end
    end
    
    directoryOut = [dirOut '/v' ID '_' methodName '_SUMMARY/'];
    [status,message,~] = mkdir(directoryOut);
    
    for fp = 1:length(fP)
        num = 0;
        ll = dir(dirOut);
        ll = {ll.name};
        
        for epoch  = 1:EEG.trials
            if strcmp(methodName,'eMI')
                if sum(strcmp(ll,['v' ID '_' methodName '_epoch' num2str(epoch) '_results.mat']))>0
                    load([dirOut 'v' ID '_' methodName '_epoch' num2str(epoch) '_results.mat'])
                    if num == 0
                        MIxtime = zeros(EEG.trials,length(fA));
                        MIpvalxtime = NaN(EEG.trials,length(fA));
                        MIStatxtime = zeros(EEG.trials,length(fA));
                        Origins     = zeros(EEG.trials,length(fA));
                    end
                    fpidx = find(fP==fP(fp));
                    if plotWithoutMask
                        MIxtime(epoch,:) = comodulogram(:,fpidx)';
                    end
                    if plotWithMask
                        MIpvalxtime(epoch,:) = pval(:,fpidx)';
                        MIStatxtime(epoch,:) = comodulogramStat(:,fpidx)';
                        Origins(epoch,:) = couplingOrigins(:,fpidx)';
                    end
                    num = 1;
                end
            else
                if sum(strcmp(ll,['Epoch' num2str(epoch)]))>0
                    load([dirOut 'Epoch' num2str(epoch) '/' 'comodulogram.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'comodulogramStat.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'pval.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'fP.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'fA.mat'])
                    if num == 0
                        MIxtime = zeros(EEG.trials,length(fA));
                        MIpvalxtime = NaN(EEG.trials,length(fA));
                        MIStatxtime = zeros(EEG.trials,length(fA));
                    end
                    fpidx = find(fP==fP(fp));
                    if plotWithoutMask
                        MIxtime(epoch,:) = comodulogram(:,fpidx)';
                    end
                    if plotWithMask
                        MIpvalxtime(epoch,:) = pval(:,fpidx)';
                        MIStatxtime(epoch,:) = comodulogramStat(:,fpidx)';
                    end
                    num = 1;
                end
            end
        end

        Time = EEG.pnts/EEG.srate*(1:1:EEG.trials);
        units = 's';
        if plotWithMask
            [~, FDRmask] = fdr(MIpvalxtime,alpha,'nonParametric');
            afterFDR = (MIxtime.*FDRmask)';
            if strcmp(methodName,'eMI')
                origins = (Origins.*FDRmask)';
                afterFDR(origins==0) = afterFDR(origins==0)*-1;
            end
        end
        if plotWithoutMask
            noFDR = MIxtime';
        end

        %% ------------------- PLOTING AND SAVING ----------------------
        if (EEG.pnts/EEG.srate*EEG.trials>60)
            Time = Time/60;
            units = 'min';
        end
        
        Titles = {};
        if plotWithMask
            Titles = [Titles 'Stat'];
        end
        if plotWithoutMask
            Titles = [Titles ''];
        end
        for i = 1:length(Titles)
            title = Titles{i};
            switch title
                case ''
                    img = noFDR;
                    cm = makeColormap(64);
                    cm = cm(ceil(size(cm,1)/2):end,:);
                    if sum(sum(img~=0))>0
                        cax = [0 max(max(img))];
                    else
                        cax = [0 1];
                    end
                case 'Stat'
                    img = afterFDR;
                    if strcmp(methodName,'eMI')
                        if sum(sum(img~=0))>0
                            cax = [-max(max(abs(img))) max(max(abs(img)))];
                        else
                            cax = [-1 1];
                        end
                        cm = makeColormap(64);
                    else
                        cm = makeColormap(64);
                        cm = cm(ceil(size(cm,1)/2):end,:);
                        if sum(sum(img~=0))>0
                            cax = [0 max(max(img))];
                        else
                           cax = [0 1]; 
                        end
                    end
                    
            end
            % rysowanie i podstawowe ustawienia osi
            f = figure('visible','off');
            set(f,'DefaultAxesFontName', fontNm)
            set(f,'DefaultAxesFontSize', fontSize)
            imagesc(Time,fA,img)
            colormap(cm);
            caxis(cax);
            set(gca, 'Ydir', 'normal')
            set(gca,'TickDir','out')
            box off
            xlabel(['Time [' units ']'])
            ylabel('Frequency for amplitude [Hz]')

            % zwezanie mapki zeby colorbar sie zmiescil
            pos = get(gca,'Position');
            set(gca,'Position',[pos(1) pos(2) pos(3)-0.05 0.98*pos(4)])

            % ustawienia colorbar
            pos = get(gca,'Position');
            axes('Position',[0.98*(pos(1)+pos(3)) pos(2) 0 pos(4)]);
            set(gca,'visible','off')
            c = colorbar();
            caxis(cax);
            cWidth = get(c,'Position');
            cWidth(3) = cWidth(3)*1.1;
            set(c,'Position',cWidth)
            set(c,'TickDir','out')
            box off
            
            if strcmp(title,'Stat')
                if strcmp(methodName,'eMI')
                    axes('Position',[1.01*(pos(1)+pos(3)) pos(2)-0.02 0.07 pos(4)*0.07]);
                    set(gca,'visible','off')
                    text(0,0,'Ambiguous','FontName',fontNm,'FontSize',fontSize*0.8,'Color',[0.3 0.3 0.3])
                    box off
                    axes('Position',[1.01*(pos(1)+pos(3)) 1.02*pos(2)+pos(4)+0.02  0.07 pos(4)*0.07]);
                    set(gca,'visible','off')
                    text(0,0,'Reliable','FontName',fontNm,'FontSize',fontSize*0.8,'Color',[1,118,202]/255)
                    box off
                end
            end

            print(f, '-dpng', '-r300', [directoryOut num2str(fP(fp)) 'Hz_' title '.png']);
            close(f);
        end
    end

end
