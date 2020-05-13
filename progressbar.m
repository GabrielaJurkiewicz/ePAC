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

function [] = progressbar(arg)

persistent buffer;
strLength  = 5;
strDotsMax = 100;

if strcmp(arg,'on')
    buffer = [];
elseif strcmp(arg,'off')
    fprintf('\n');
    clear buffer
elseif isnumeric(arg)
    arg = floor(arg);
    percentageOut = [num2str(arg) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strLength-length(percentageOut)-1)];
    nDots = floor(arg/100*strDotsMax);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMax-nDots) ']'];
    strOut = [percentageOut dotOut];
    fprintf([buffer strOut]);
    buffer = repmat('\b',1,length(strOut)-1);
end

end