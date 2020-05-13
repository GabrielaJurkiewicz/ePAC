%    Phase-amplitude coupling detection MATLAB plugin
%
%    Copyright (C) 2012 Jaroslaw Zygierewicz & Maciej Kaminski
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

%    Jaroslaw Zygierewicz <Jaroslaw.Zygierewicz@fuw.edu.pl>
%    Maciej Kaminski <Maciek.Kaminski@fuw.edu.pl>

function output=tf_morlet(M, varargin)

% Morlet wavelet
% M - length of the wavelet
% w - if present, Morlet wavelet paramter(number of oscillations per
% wavelet). Default value w = 5.
% s - if present, scaling factor. Default value s = 1.
% complete - if present, whether to use the complete or the standard version. Default value complete = True.
% The complete version of the Morlet wavelet, with a correction term
% improves admissibility. For w greater than 5, the correction term is negligible.
% The admissibility is a condition for a successful inverse transform and
% implies that a wavelet must integrate to zero.

w = 5;
s = 1;
complete = true;
switch nargin
    case 0
        disp('tf_morlet - not enough arguments given');
        return
    case 1

    case 2
        w=varargin{1};
    case 3
        w=varargin{1};
        s=varargin{2};
    otherwise
        w=varargin{1};
        s=varargin{2};
        complete=varargin{3};
end

x = linspace(-s*2*pi, s*2*pi, M);
output = exp(i*w*x);

if complete
    output = output-exp(-0.5*w^2);
end

output = output.*exp(-0.5*(x.^2)).*pi^(-0.25);
