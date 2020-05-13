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

function [] = plot_maps_together(x,y,fSpectrum,mapa,signal,meanTheta,CouplingOrig,PC,PhaseBins,SA,AS,comod,w,Fs,savename,cut,plotWithMask)

    AS = AS/sum(AS)*100;
    SA = SA/sum(SA)*100;
    
    if plotWithMask
        phase = angle(hilbert(meanTheta));
        mask  = zeros(size(mapa));
        for fa = 1:length(y)
            IDX = zeros(size(phase));
            PH = find(PC(fa,:)>0);
            if ~isempty(PH)
                for ph = 1:length(PH)
                    idxPH = ((phase>PhaseBins(PH(ph))) & (phase<=PhaseBins(PH(ph)+1)));
                    IDX = IDX | idxPH;
                end
            end
            mask(fa,:) = IDX;
        end
        [fmx,fwd,fstart,fstop,fspmx,spmx] = findMaxInComodulogram(comod,SA,AS,fSpectrum,y,w,Fs);
    else
        fmx = [];
    end
   
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
    labelSize = 10;
    f = figure('visible','off');
    set(f,'DefaultAxesFontName', fontNm)
    set(f,'DefaultAxesFontSize', labelSize)
    
    subplot(5,1,[1,2,3,4]);
    imagesc(x,y,log(mapa))
    cm = makeColormap(64);
    colormap(cm(ceil(size(cm,1)/2):end,:))
    if sum(sum(mapa~=0))>0
        caxis([min(min(log(mapa))) max(max(log(mapa)))])
    else
        caxis([0 1])
    end
    title('Averaged time-frequency map','FontSize',labelSize,'FontName',fontNm)
    if plotWithMask
        if sum(sum(mask))>0
            hold on;
            contour(x,y,mapa.*mask,[eps eps],'.','Color',[0 0 0]/255,'LineWidth',1.5);
        end
    end
    set(gca,'position',[0.23,0.35,0.55,0.60])
    set(gca,'Ydir','normal')
    set(gca,'XTick',[])
    set(gca,'TickDir','out')
    box off
    ylabel('Frequency [Hz]','FontSize',labelSize,'FontName',fontNm)

    c = colorbar();
    set(c,'position',[0.1,0.35,0.03,0.60])
    set(c,'TickDir','out')
    if sum(sum(mapa~=0))>0
        caxis([min(min(log(mapa))) max(max(log(mapa)))])
    else
        caxis([0 1])
    end
    box off
    ylabel(c,'log(power)','FontSize',labelSize,'FontName',fontNm)
    
    axx = axes();
    set(axx,'position',[0.79,0.35,0.15,0.60])
    if (sum(AS>0)>0)&&(sum(SA>0)>0)
        xmin = min([AS; SA]);
        xmax = max([AS; SA]);
        xlim([xmin xmax])
        ylim([min(fSpectrum) max(fSpectrum)])
        if ~isempty(fmx)
            for fm = 1:length(fmx)
                if sum(CouplingOrig((y>=fstart(fm))&(y<=fstop(fm))))>0
                    col = [1,184,202]/255;
                else
                    col = [0 0 0];
                end
                fSp = fSpectrum((fSpectrum>=fmx(fm)-fwd(fm)/2)&(fSpectrum<=fmx(fm)+fwd(fm)/2));
                fill([xmin xmin xmax xmax],[min(fSp) max(fSp) max(fSp) min(fSp)],col,'FaceAlpha',0.2,'EdgeColor','none')
                hold on
            end
        end
    else
        xmin = 0;
        xmax = 1;
    end
    plot(AS,fSpectrum,'Color',[0 0 0],'LineWidth',1.5)
    hold on
    plot(SA,fSpectrum,'Color',[1,184,202]/255,'LineWidth',1.5)
    ylim([min(fSpectrum) max(fSpectrum)])
    xlim([xmin xmax])
    if ~isempty(fmx)
        for fm = 1:length(fmx)
            hold on
            plot(spmx(fm),fspmx(fm),'o','Color',[1 108 202]/255,'LineWidth',1.5,'MarkerSize',4)
        end
    end
    xlabel('power [%]','FontSize',labelSize,'FontName',fontNm)
    set(gca,'YTickLabel',[])
    set(gca,'YTick',[])
    set(gca,'TickDir','out')
    box off
    
    axes('position',[0.93,0.35,0.05,0.60]);
    set(gca,'visible','off')
    xlim([0 1])
    ylim([0 1])
    text(0.5,1,'Spectrum of averaged signal      and average spetrum   ','FontSize',0.8*labelSize,'FontName',fontNm,'FontWeight','bold','Rotation',270)
    hold on
    plot([0.5 0.5],[0.45 0.48],'Color',[1,184,202]/255,'LineWidth',2)
    hold on
    plot([0.5 0.5],[0.03 0],'Color',[0,0,0]/255,'LineWidth',2)
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    box off
 
    subplot(5,1,5);
    plot(x,signal,'Color',[1,184,202]/255,'LineWidth',2)
    hold on
    plot(x,meanTheta,'Color',[0,0,0]/255,'LineWidth',1)
    set(gca,'position',[0.23,0.11,0.55,0.17])
    set(gca,'XTick',round(linspace(-cut,cut,7),2,'significant'))
    set(gca,'TickDir','out')
    box off
    xlim([-cut, cut])
    xlabel('Time [s]','FontSize',labelSize,'FontName',fontNm)
    ylabel('Amplitude','FontSize',labelSize,'FontName',fontNm)
    
    axes('position',[0.23,0.29,0.55,0.05]);
    set(gca,'visible','off')
    xlim([0 1])
    ylim([0 1])
    text(0,0.5,'Averaged high-freq signal     and low-freq oscillation    ','FontSize',labelSize,'FontName',fontNm,'FontWeight','bold')
    hold on
    plot([0.48 0.51],[0.5 0.5],'Color',[1,184,202]/255,'LineWidth',2)
    hold on
    plot([0.97 1],[0.5 0.5],'Color',[0,0,0]/255,'LineWidth',2)
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    box off
    
    print(f, '-dpng', '-r300', [savename '.png']);
    close(f);
    
end