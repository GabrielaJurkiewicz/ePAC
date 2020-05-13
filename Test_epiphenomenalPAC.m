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
fPstart    = 2;
fPend      = 30;
fPstep     = 1;
fAstart    = 30;
fAend      = 100;
fAstep     = 5;
Nboot      = 200;
artifacts  = true;
plotWithMask    = 1;
plotWithoutMask = 0;


%% EXAMPLE #1a OF EPIPHENOMENAL COUPLING: Periodic trains of Gaussian functions
%    This example display the results of analysis of signal containing epiphenomenal phase-amplitude coupling. This signal is a replication of the one 
%    proposed by Gerber et al. (2016) in "Nonsinusoidal activity can produce cross-frequency coupling in cortical signals in the absence of functional 
%    interaction between neural sources". It is used to demonstrate how a semi-periodic occurrence of sharp waveforms can drive epiphenomenal coupling. 
%    In this case, the slow-rhythm is produced by a train of Gaussian-shaped spikes placed at intervals drawn from a uniform distribution of 100+-20 ms, 
%    which corresponds to spectral peak at 10 Hz. This signal is added to an EEG signal that does not contain any PAC. The Gaussian spikes amplitude is 
%    set to 2 standard deviations of the background EEG signal. After running the program the output images and files will be produced. In the picture 
%    v1_eMI_comodulogramStat.png you can see that, as expected, the algorithm detects coupling labeled as ambiguous for frequency for phase ~10 Hz and 
%    for wide range of frequencies for amplitude. In the picture  v1_eMI_results10Hz.png we present the auxiliary plots for fP=10 Hz specifically and 
%    we can observe that coupling occurs where the average signal undergoes abrupt changes in amplitude. Also, we can note that the averaged signal 
%    (turquoise line, bottom plot) does not contain a sinusoidal low-frequency component as extracted through filtration (black line, bottom plot), 
%    and the time-frequency maps display a wide range of coupled frequencies. The average spectrum (black line) does not contain any prominent peaks, 
%    whereas the spectrum of the averaged signal (turquoise line) displays characteristic harmonic structure with a maximum not congruent with the 
%    maximum in comodulogram. All of those are an indicator of epiphenomenal origins. The results also include the image v1_eMI_phaseHistogramStat.png 
%    where we can observe the coupled phase. The reference dPAC (v1_dPAC_comodulogramStat.png) and MI (v1_MI_comodulogramStat.png) methods 
%    provide only comodulograms which display results consistent with eMI method.

EEG = pop_loadset('filename', 'Test_epiphenomenalPAC_1a.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_epiphenomenalPAC_1a/';

% running Extended Modulation Index Analysis
fPband     = 2;
w          = 5;
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

        
%% EXAMPLE #1b OF EPIPHENOMENAL COUPLING: Periodic trains of Gaussian functions
%    This signal is analogous to the example #1a except the Gaussian spikes amplitude is set to 5 standard deviations of the background EEG signal. 
%    Here the results are also analogous to the #1a example except all of the indicators of epiphenomenal origins are much more pronounced. The 
%    reference dPAC (v1_dPAC_comodulogramStat.png) and MI (v1_MI_comodulogramStat.png) methods provide only comodulograms which display results 
%    consistent with eMI method.

EEG = pop_loadset('filename', 'Test_epiphenomenalPAC_1b.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_epiphenomenalPAC_1b/';

% running Extended Modulation Index Analysis
fPband     = 2;
w          = 5;
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
        
    
%% EXAMPLE #2a OF EPIPHENOMENAL COUPLING: Non-periodic trains of Gaussian functions
%    This example also display the results of analysis of signal containing epiphenomenal phase-amplitude coupling. This signal is used to demonstrate 
%    that periodicity is not strictly necessary for PAC as long as the sharp waveforms are spaced sufficiently apart for a given frequency-for-phase 
%    such that the peaks of this slow frequency can align with the sharp waves (Gerber et al. 2016, "Nonsinusoidal activity can produce cross-frequency 
%    coupling in cortical signals in the absence of functional interaction between neural sources"). The signal was obtained in the following way. 
%    A series of 100 random time values are obtained from a uniform distribution ranging from 0 to 10 s with 0.001 s step. A Gaussian spike is placed 
%    at each event, and the resulting signal is superimposed with a 10 s background EEG signal that does not contain any PAC. The Gaussian spike amplitude 
%    is set to 2 standard deviations of the background EEG signal. After running the program the output images and files will be produced. In the picture 
%    v1_eMI_comodulogramStat.png you can see that the signal did not produce any statistically significant PAC. The reference dPAC (v1_dPAC_comodulogramStat.png) 
%    and MI (v1_MI_comodulogramStat.png) methods provide only comodulograms which display results consistent with eMI method.

EEG = pop_loadset('filename', 'Test_epiphenomenalPAC_2a.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_epiphenomenalPAC_2a/';

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
        
        
%% EXAMPLE #2b OF EPIPHENOMENAL COUPLING: Non-periodic trains of Gaussian functions
%    This signal is analogous to the example #2a except the Gaussian spikes amplitude is set to 5 standard deviations of the background EEG signal. 
%    After running the program the output images and files will be produced. In the picture v1_eMI_comodulogramStat.png you can see that most of the 
%    coupling regions are correctly labelled as ambiguous. However, there is one region with falsely assigned Reliable label, which should be a reminder 
%    that eMI also has some limitations. The automatic label assignment sometimes may be wrong, which is why we strongly suggest double-check it by visual
%    inspection of auxiliary plots. In the picture v1_eMI_results2Hz.png we present the auxiliary plots for fP=2 Hz specifically and we can observe that 
%    the coupling occurs where the average signal undergoes abrupt changes in amplitude. Also, we can note that the averaged signal does not contain a 
%    sinusoidal low-frequency component as extracted through filtration, and the time-frequency maps display a wide range of coupled frequencies. 
%    The average spectrum does not contain any prominent peaks, whereas the spectrum of the averaged signal display characteristic harmonic structure. 
%    All of those features indicate the epiphenomenal origins of the coupling. The results also include the image v1_eMI_phaseHistogramStat.png 
%    where we can observe the coupled phase. The reference dPAC (v1_dPAC_comodulogramStat.png) and MI (v1_MI_comodulogramStat.png) methods 
%    provide only comodulograms which display results consistent with eMI method. 

EEG = pop_loadset('filename', 'Test_epiphenomenalPAC_2b.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);
dirOut = './Results_epiphenomenalPAC_2b/';

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
