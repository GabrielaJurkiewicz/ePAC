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

function plotFigurePhase(fP,fA,PHASE,PhaseComod,CouplingOrigins,clabel,dirOut,savename)

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
    
    PhaseComod(PhaseComod<0) = 0;
    
    if ~isempty(CouplingOrigins)
        if sum(sum(sum(PhaseComod>0)))>0
            cax = [-max(max(max(PhaseComod))) max(max(max(PhaseComod)))];
        else
            cax = [-1 1];
        end
        cm = makeColormap(64);
        for i = 1:size(PhaseComod,1)
            for j = 1:size(PhaseComod,3)
                if CouplingOrigins(i,j)==0
                    PhaseComod(i,:,j) = PhaseComod(i,:,j)*-1;
                end
            end
        end
    else
        cm = makeColormap(64);
        cm = cm(ceil(size(cm,1)/2):end,:);
        if sum(sum(sum(PhaseComod>0)))>0
            cax = [0 max(max(max(PhaseComod)))];
        else
            cax = [0 1];
        end
    end
    
    PhaseComod = reshape(PhaseComod,[length(fA),length(fP)*length(PHASE)]);

    % rysowanie i podstawowe ustawienia osi
    f = figure('visible','off');
    set(f,'DefaultAxesFontName', fontNm)
    set(f,'DefaultAxesFontSize', fontSize)
    imagesc(1:1:length(fP)*length(PHASE),fA,PhaseComod)
    colormap(cm);
    caxis(cax);
    xlim([0.5 (length(fP)*length(PHASE))+1])
    set(gca, 'Ydir', 'normal')
    set(gca,'TickDir','out')
    
    % zwezanie mapki zeby colorbar sie zmiescil
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) pos(3)-0.05 0.98*pos(4)])
    
    % rysowanie pionowych bialych linii rozdzielacjacych fP
    for l = 2:length(fP)
        hold on 
        line([(l-1)*length(PHASE)+0.5 (l-1)*length(PHASE)+0.5],get(gca,'YLim'),'LineWidth',1.5,'LineStyle',':','Color',[227 219 219]/255)%[134 118 220]/255
    end
    
    % ustawienia osi X i Y mapki; X faza wyswietlana na gorze mapy
    phase_ticks = repmat(PHASE',length(fP),1)';
    phase_labels = {};
    idxP = find((phase_ticks == -pi) | (phase_ticks == 0));
    for pt = 1:length(idxP)
        if phase_ticks(idxP(pt))==-pi
            if pt==1
                phase_labels = [phase_labels '-\pi'];
            else
                phase_labels = [phase_labels '\pi'];
            end
        else
            phase_labels = [phase_labels '0'];
        end
    end
    idxP = [idxP length(fP)*length(PHASE)+1];
    phase_labels = [phase_labels '\pi'];
    set(gca,'XTick',idxP-0.5)
    set(gca,'XTickLabel',phase_labels)
    set(gca,'xaxisLocation','top')
    box off
    xlabel('Phase')
    ylabel('Frequency for amplitude [Hz]')
    
    % ustawienia colorbar
    pos = get(gca,'Position');
    axes('Position',[0.98*(pos(1)+pos(3)) pos(2) 0 pos(4)]);
    set(gca,'visible','off')
    c = colorbar();
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
    
    % dodatkowa os na dole mapy zawierajaca podpisy czestosci fP
    fPlabel = {};
    for p = 1:length(fP)
        fPlabel = [fPlabel num2str(fP(p))];
    end
    axes('Position',[pos(1) pos(2) pos(3) 0]);
    set(gca,'XTick',((length(PHASE)/(length(PHASE)*length(fP)):length(PHASE)/(length(PHASE)*length(fP)):1)-length(PHASE)/(length(PHASE)*length(fP)*2)))
    set(gca,'XTickLabel',fPlabel)
    set(gca,'TickDir','out')
    set(gca,'xaxisLocation','bottom')
    set(gca,'TickLength',[0 0])
    xlabel('Frequency for phase [Hz]')
    box off
    
    print(f, '-dpng', '-r300', [dirOut savename '.png']);
    close(f);
    
end