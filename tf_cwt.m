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

function [P,freqs,t]=tf_cwt(x,MinF,MaxF,Fs,varargin)

% Continuous wavelet transform
% x - signal
% MinF, MaxF - frequency range [Hz]
% Fs  - sampling frequency
% w - if present, Morlet wavelet paramter(number of oscillations per
% wavelet). Default value w = 7.
% df - if present, frequency step. Default value df = 7.

w = 7;
df = 1;
% doplot = true;
switch nargin
    case 0:3
        disp('tf_cwt - not enough input arguments');
        return
    case 4

    case 5
        w=varargin{1};
    case 2
        w=varargin{1};
        df=varargin{2};
    otherwise
        w=varargin{1};
        df=varargin{2};
%         doplot=varargin{3};
end

T = length(x)/Fs;
M = length(x);
t = 0:1/Fs:T; t(end) = [];
freqs = MinF:df:MaxF; freqs(end) = [];
P = zeros(length(freqs),M);
X = fft(x);

for i=1:length(freqs)
    f = freqs(i);
    s = T*f/(2*w);
    psi = fft(tf_morlet(M,w,s,true));
    psi = psi./sqrt(sum(abs(psi).^2));
    tmp = fftshift(ifft(X.*psi));
    P(i,:) = abs(tmp).^2;
end

% if doplot
%     pcolor(linspace(0,T,size(P,2)),linspace(MinF,MaxF,size(P,1)),P); shading interp;
% end