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

%% SETTING COMMON PARAMETERS
eeglab

chan_fP    = 1;
chan_fA    = 1;
epochs     = [1];
type       = 'continuous';
fPstart    = 4;
fPend      = 8;
fPstep     = 1;
fAstart    = 20;
fAend      = 120;
fAstep     = 5;
Nboot      = 200;
artifacts  = true;
plotWithMask    = 1;
plotWithoutMask = 0;


%% EXAMPLE #1 OF PROPER COUPLING: coupled oscillatory bursts model
%    This example display the results of analysis of artificial signal with proper phase-amplitude coupling. It consists of 6 Hz sinus with 
%    superimposed bursts of 77 Hz oscillation located just before the top of each sinus cycle with added white noise. This synthetic signal 
%    represents the PAC where phase of 6 Hz oscillation modulates the amplitude of 77 Hz osillation. After running the program the output images 
%    and files will be produced in the indicated folder. In the picture v1_eMI_comodulogramStat.png you can see that the algorithm correctly detects 
%    significant coupling labeled as reliable for frequency for phase = 6 Hz and for the range of frequencies for amplitude around 77 Hz. In the picture 
%    v1_eMI_results6Hz.png we present the auxiliary plots for fP=6 Hz specifically and you can observe that the augmented and outlined regions in the 
%    time-frequency map indicate the phase of low-frequency oscillation (bottom plot, black line) which is just before the top. That is in line with 
%    the PAC that was simulated. Additionally you can obsereve the prominent peak around 77Hz in average spectrum (black line) and spectrum of averaged signal
%    (turquoise line). The results also include the image v1_eMI_phaseHistogramStat.png where we can observe that the coupling occurs for the range 
%    of phases consistent with the simulation parameters. The reference dPAC (v1_dPAC_comodulogramStat.png) and MI (v1_MI_comodulogramStat.png) methods 
%    provide only comodulograms which display results consistent with eMI method.

EEG = pop_loadset('filename', 'Test_properPAC_1.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_properPAC_1/';

% running Extended Modulation Index Analysis
fPband     = 1;
w          = 5;
nbBins     = 18;
peMI       = 95;
pPhaseCom  = 95;
eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,...
             pPhaseCom,artifacts,plotWithMask,plotWithoutMask);

% running Modulation Index Analysis
fPfiltband = 1;
fAfiltband = fPend*2;
pMI        = 95;
MI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,nbBins,Nboot,pMI,...
            artifacts,plotWithMask,plotWithoutMask);

% running direct PAC Estmate Analysis
fPfiltband = 1;
fAfiltband = fPend*2;
pdPAC      = 95;
dPAC_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,Nboot,pdPAC,...
            artifacts,plotWithMask,plotWithoutMask);
        
        
%% EXAMPLE #2 OF PROPER COUPLING: amplitude modulation model
%    This exmaple display the results of analysis of another artificial signal with proper phase-amplitude coupling. It consists of 6 Hz sinus, 
%    77 Hz sinus with amplitude altering according to the phase of the low-frequency sine (maximal amplitude for the phase 5*pi/3). This synthetic signal 
%    represents the PAC where phase of 6 Hz oscillation modulates the amplitude of 77 Hz osillation. After running the program the output images and files 
%    will be produced in the indicated folder. In the picture v1_eMI_comodulogramStat.png you can see that the algorithm correctly detects significant 
%    coupling labeled as reliable for frequency for phase = 6 Hz and for the range of frequencies for amplitude around 77 Hz. In the picture 
%    v1_eMI_results6Hz.png we present the auxiliary plots for fP=6 Hz specifically and you can observe that the augmented and outlined regions in the 
%    time-frequency map indicate the phase of low-frequency oscillation (bottom plot, black line) which is just before the top. That is in line with the 
%    PAC that was simulated. Additionally you can obsereve the prominent peak around 77Hz in average spectrum (black line) and spectrum of averaged signal 
%    (turquoise line). The results also include the image v1_eMI_phaseHistogramStat.png where we can observe that the coupling occurs for the range of 
%    phases consistent with the simulation parameters. The reference dPAC (v1_dPAC_comodulogramStat.png) and MI (v1_MI_comodulogramStat.png) methods 
%    provide only comodulograms which display results consistent with eMI method.

EEG = pop_loadset('filename', 'Test_properPAC_2.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_properPAC_2/';

% running Extended Modulation Index Analysis
fPband     = 1;
w          = 5;
nbBins     = 18;
peMI       = 95;
pPhaseCom  = 95;
eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,...
             pPhaseCom,artifacts,plotWithMask,plotWithoutMask);

% running Modulation Index Analysis
fPfiltband = 1;
fAfiltband = fPend*2;
pMI        = 95;
MI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,nbBins,Nboot,pMI,...
            artifacts,plotWithMask,plotWithoutMask);

% running direct PAC Estmate Analysis
fPfiltband = 1;
fAfiltband = fPend*2;
pdPAC      = 95;
dPAC_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,Nboot,pdPAC,...
            artifacts,plotWithMask,plotWithoutMask);
