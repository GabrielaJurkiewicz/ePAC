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

function [] = plotPolarHist(fP,fA,phaseBins,comodulogram,phaseComodulogram,CouplingOrigins,nbBins,dirOut,savename)

    fontSize = 10;
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
    
    dp = (2*pi/nbBins)/2;
    [Phase,Colors,Orig] = preparePhaseHist(fP,fA,phaseBins+dp,comodulogram,phaseComodulogram,CouplingOrigins);
    
    if ~isempty(Phase)
        f = figure('visible','off');
        set(f,'DefaultAxesFontName', fontNm)
        set(f,'DefaultAxesFontSize', fontSize)
        v = version('-release');
        vYear = str2num(v(1:4));
    
        axes('Position',[0.1 0.3 0.8 0.6])
        if vYear>2016 || strcmp(v,'2016b')
            for r = 1:length(Phase)
                pp = polarhistogram(Phase{end-r+1},[phaseBins pi],'Normalization','probability','DisplayStyle','stairs');
                set(pp,'EdgeColor',Colors{end-r+1})
                switch Orig{end-r+1}
                    case 'reliable'
                        set(pp,'LineWidth',3)
                    case 'ambiguous'
                        set(pp,'LineWidth',4)
                end
                hold on
            end
            set(pp.Parent,'ThetaAxisUnits','radians')
            set(pp.Parent,'ThetaZeroLocation','right')
            set(pp.Parent,'Box','off')
            set(pp.Parent,'FontName',fontNm)
            set(pp.Parent,'FontSize',fontSize)
            set(pp.Parent,'TickDir','out')
            
            axes('Position',[0.3 0.1 0.4 0.1])
            tim = 0:1/500:2.5;
            plot(tim,cos(2*pi*1*tim+pi/2),'Color',[0 0 0])
            xlim([0 1.5])
            set(gca,'TickDir','out');box off;set(gca,'YTick',[]);set(gca,'TickLength',[0.04 0.04])
            set(gca,'XTick',0:0.25:1.5)
            set(gca,'XTickLabel',{'','-\pi','-\pi/2','0','\pi/2','\pi',''})
            set(gca,'YColor','none','XColor',[0.5 0.5 0.5])
            set(gca,'FontSize',fontSize)
            set(gca,'FontName',fontNm)
            
        else
            for r = 1:length(Phase)
                [tt, rr] = rose(Phase{end-r+1},phaseBins+dp);
                rr = rr./numel(Phase{end-r+1});
                pp = polar(tt, rr);
                set(pp,'Color',Colors{end-r+1})
                switch Orig{end-r+1}
                    case 'reliable'
                        set(pp,'LineWidth',3)
                    case 'ambiguous'
                        set(pp,'LineWidth',4)
                end
                hold on
            end
            set(pp.Parent,'Box','off')
            set(pp.Parent,'FontName',fontNm)
            set(pp.Parent,'FontSize',fontSize)
            set(pp.Parent,'TickDir','out')
            
            axes('Position',[0.3 0.1 0.4 0.1])
            tim = 0:1/500:2.5;
            plot(tim,cos(2*pi*1*tim+pi/2),'Color',[0 0 0])
            xlim([0 1.5])
            set(gca,'TickDir','out');box off;set(gca,'YTick',[]);set(gca,'TickLength',[0.04 0.04])
            set(gca,'XTick',0:0.25:1.5)
            set(gca,'XTickLabel',{'','-180','-90','0','90','180',''})
            set(gca,'YColor','none','XColor',[0.5 0.5 0.5])
            set(gca,'FontSize',fontSize)
            set(gca,'FontName',fontNm)
        end

        print(f, '-dpng', '-r300', [dirOut savename '.png']);
        close(f);
    end
    
end