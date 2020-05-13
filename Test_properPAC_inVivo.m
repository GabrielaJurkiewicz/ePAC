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

%% EXAMPLE #3 OF PROPER COUPLING: In vivo recording from rat after ketamine injection
%    This example display the results of analysis of the in vivo Local Field Potentials recording from rat after ketamine injection. This data come 
%    from the study "The olfactory bulb is a source of high-frequency oscillations 130-180 Hz associated with a subanesthetic dose of ketamine 
%    in rodents", Hunt MJ, Adams NE, Œredniawa W, Wójcik DK, Simon A, Kasicki S, Whittington MA, Neuropsychopharmacology 44(2):435–442, 
%    DOI 10.1038/s41386-018-0173-y. It is a section (1560-1580 sec) of recording from the right side of olfactory bulb of the RAT86, 6 minutes 
%    after the injection of ketamine, which is expected to induce the High Frequency Oscillations (around 150 Hz) coupled with theta rhytm.   
%    After running the program the output images and files will be produced in the indicated folder. In the picture v1_eMI_comodulogramStat.png 
%    you can see that the algorithm detects expectd coupling labeled as reliable for frequency for phase ~ 7 Hz and for the range of frequencies 
%    for amplitude around 150 Hz. In the picture v1_eMI_results7Hz.png we present the auxiliary plots for fP=7 Hz specifically and you can observe 
%    the augmented and outlined regions in the time-frequency map which indicate coupling with the phase of the low-frequency oscillation 
%    (top of the theta cycle). Additionally you can obsereve the prominent peak around 150 Hz in average spectrum (black line) and spectrum of averaged signal
%    (turquoise line). The results also include the image v1_eMI_phaseHistogramStat.png where we can observe that the coupling occurs for the phase ~ 0 pi
%    (top of the thata cycle). The reference dPAC (v1_dPAC_comodulogramStat.png) and MI (v1_MI_comodulogramStat.png) methods provide only comodulograms 
%    which display results consistent with eMI method.

eeglab

chan_fP    = 1;
chan_fA    = 1;
epochs     = [1];
type       = 'continuous';
fPstart    = 2;
fPend      = 10;
fPstep     = 1;
fAstart    = 40;
fAend      = 190;
fAstep     = 5;
Nboot      = 200;
artifacts  = true;
plotWithMask    = 1;
plotWithoutMask = 0;

EEG = pop_loadset('filename', 'Test_properPAC_inVivo.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_properPAC_inVivo/';

% running Extended Modulation Index Analysis
fPband     = 2;
w          = 7;
nbBins     = 18;
peMI       = 95;
pPhaseCom  = 95;
eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,...
             pPhaseCom,artifacts,plotWithMask,plotWithoutMask);

% running Modulation Index Analysis
fPfiltband = 2;
fAfiltband = fPend*2;
pMI        = 95;
MI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,nbBins,Nboot,pMI,...
            artifacts,plotWithMask,plotWithoutMask);

% running direct PAC Estmate Analysis
fPfiltband = 2;
fAfiltband = fPend*2;
pdPAC      = 95;
dPAC_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,Nboot,pdPAC,...
            artifacts,plotWithMask,plotWithoutMask);
