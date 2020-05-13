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

function [Phase,Colors,Orig] = preparePhaseHist(fP,fA,phaseBins,comodulogram,phaseComodulogram,CouplingOrigins)

    Phase = {};
    Orig = {};
    [~,~,Colors,L_rel,L_ambig] = findBoundaries(fP,fA,comodulogram,CouplingOrigins);
    nums_rel = unique(L_rel(L_rel~=0));
    nums_ambig = unique(L_ambig(L_ambig~=0));
    for r = 1:length(nums_rel) 
        phase = [];
        loc = (L_rel==nums_rel(r));
        for fa = 1:length(fA)
            for fp = 1:length(fP)
                if loc(fa,fp)==1
                    phase = [phase phaseBins(squeeze(phaseComodulogram(fa,:,fp))>0)];
                end
            end
        end
        Phase = [Phase, phase]; 
        Orig = [Orig 'reliable'];
    end
    for a = 1:length(nums_ambig)
        phase = [];
        loc = (L_ambig==nums_ambig(a));
        for fa = 1:length(fA)
            for fp = 1:length(fP)
                if loc(fa,fp)==1
                    phase = [phase phaseBins(squeeze(phaseComodulogram(fa,:,fp))'>0)];
                end
            end
        end
        Phase = [Phase, phase]; 
        Orig = [Orig 'ambiguous'];
    end

end