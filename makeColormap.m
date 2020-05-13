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

function [cm] = makeColormap(M)

    if rem(M,2)==0
        m = M/2;
    else
        m = (M-1)/2;
    end
    bottom = [242,253,255]/255;
    middle1 = [210 246 249]/255;
    middle2 = [46 200 216]/255;
    top = [1 118 202]/255;
    
    new = [bottom; middle1; middle2; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    for i=1:3
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    cm1 = newmap;
    
    bottom = [0.20 0.20 0.20];
    middle = [0.6 0.6 0.6];
    top = [0.96,0.96,0.96];
    
    new = [bottom; middle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    for i=1:3
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    cm2 = newmap;
    
    cm = [ [cm2(:,1);1;cm1(:,1)],[cm2(:,2);1;cm1(:,2)],[cm2(:,3);1;cm1(:,3)]];

end