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

function plotFigure(fP,fA,comodulogram,CouplingOrigins,clabel,dirOut,savename)

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
    
    comodulogram(comodulogram<0) = 0;
    
    if ~isempty(CouplingOrigins)
        if sum(sum(comodulogram>0))>0
            cax = [-max(max(comodulogram)) max(max(comodulogram))];
        else
            cax = [-1 1];
        end
        comodulogram(CouplingOrigins==0) = comodulogram(CouplingOrigins==0)*-1;
        cm = makeColormap(64);
    else
        cm = makeColormap(64);
        cm = cm(ceil(size(cm,1)/2):end,:);
        if sum(sum(comodulogram>0))>0
            cax = [0 max(max(comodulogram))];
        else
            cax = [0 1];
        end
    end
    
    f = figure('visible','off');
    set(f,'DefaultAxesFontName', fontNm)
    set(f,'DefaultAxesFontSize', fontSize)
    imagesc(fP,fA,comodulogram)
    colormap(cm);
    caxis(cax);
    if ~isempty(CouplingOrigins)
        [X,Y,Colors,~,~] = findBoundaries(fP,fA,comodulogram,CouplingOrigins);
        for b = 1:length(X)
            hold on
            plot(X{end-b+1}, Y{end-b+1}, 'LineWidth', 4, 'Color', Colors{end-b+1})
        end
    end
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) pos(3)-0.05 0.98*pos(4)])
    set(gca,'Ydir','normal')
    set(gca,'TickDir','out')
    set(gca,'XTick',fP)
    set(gca,'XTickLabel',fP)
    ylabel('Frequency for amplitude [Hz]')
    xlabel('Frequency for phase [Hz]')
    box off

    pos = get(gca,'Position');
    axes('Position',[0.98*(pos(1)+pos(3)) pos(2) 0 pos(4)]);
    set(gca,'visible','off')
    c = colorbar;
    caxis(cax);
    ylabel(c,clabel)
    cWidth = get(c,'Position');
    cWidth(3) = cWidth(3)*1.1;
    set(c,'Position',cWidth)
    set(c,'TickDir','out')
    box off
    
    if ~isempty(CouplingOrigins)
        axes('Position',[1.01*(pos(1)+pos(3)) pos(2)-0.02 0.07 pos(4)*0.07]);
        set(gca,'visible','off')
        text(0,0,'Ambiguous','FontName',fontNm,'FontSize',fontSize*0.8,'Color',[0.3 0.3 0.3])
        box off
        axes('Position',[1.01*(pos(1)+pos(3)) 1.02*pos(2)+pos(4)+0.02  0.07 pos(4)*0.07]);
        set(gca,'visible','off')
        text(0,0,'Reliable','FontName',fontNm,'FontSize',fontSize*0.8,'Color',[1,118,202]/255)
        box off
    end
    
    print(f, '-dpng', '-r300', [dirOut savename '.png']);
    close(f);

end