%    extended Phase-amplitude coupling detection EEGLAB toolbox
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


function eegplugin_ePAC_calc(fig, try_strings, catch_strings)

if nargin < 3
    throw(MException('PhaseAmplitudeCoupling:eegplugin_ePAC_calc','Function eegplugin_ePAC_calc requires 3 arguments, only %s given.',len(nargin)));
end

Toolmenu = findobj(fig, 'tag', 'tools');
PACmenu = uimenu( Toolmenu, 'label', 'extended Phase-Amplitude Coupling');

com1 = ['[LASTCOM] = pop_eMI_calc(EEG);'];
com2 = ['[LASTCOM] = pop_MI_calc(EEG);'];
com3 = ['[LASTCOM] = pop_dPAC_calc(EEG);'];

uimenu( PACmenu, 'label', 'extended Modulation Index', 'callback', com1);
uimenu( PACmenu, 'label', 'Modulation Index', 'callback', com2);
uimenu( PACmenu, 'label', 'direct PAC estimate', 'callback', com3);


end
