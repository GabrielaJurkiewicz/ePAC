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
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    Gabriela Jurkiewicz <gabriela.j.jurkiewicz@gmail.com>

% -------------------------------------------------------------------------------------------------------------------------
%     eMI:
%     Calculates the Extended Modulation Index - measure of
%     phase-amplitude coupling for a given data.
% 
%     Usage:
%       >> pop_eMI_calc(EEG);          % pop_up window
%
%       >> eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbCycles,...
%                   nbBins,w,Nboot,peMI,pPhaseCom,artifacts,plotWithMask,plotWithoutMask);
% 
%     Input:
%       - Channel number for phase signal (chan_fP): 
%         number of the EEG data channel which contains the low-frequency oscillation (signal for phase)
%         that provides the phase which is the modulator of the amplitude of the high-frequency 
%         oscillation (signal for amplitude)
% 
%       - Channel number for amplitude signal (chan_fA):
%         number of the EEG data channel that contains the high-frequency oscillation (signal for amplitude)
%         which amplitude is modulated by the phase of the low-frequency oscillation (signal for phase)
%
%       - Epoch numbers (epochs):
%         numbers of epochs that will be analysed 
% 
%       - Frequency for phase start (fPstart):
%         the beginning of the examined range of frequencies for phase (frequencies of low-frequency oscillation 
%         with modulating phase) 
% 
%       - Frequency for phase stop (fPend):
%         the end of the examined range of frequencies for phase (frequencies of low-frequency oscillation  
%         with modulating phase)
% 
%       - Frequency for phase step (fPstep):
%         the step between the beginning and end of examined range of frequencies for phase. 
% 
%       - Filtration bandwidth for phase frequenies (fPband):
%         Each frequency for phase f_P is in fact a bin ranging from f_P-fPband/2 to f_P+fPband/2
%
%       - Frequency for amplitude start (fAstart):
%         the beginning of the examined range of frequencies for amplitude (frequencies of high-frequency oscillation
%         with modulated amplitude)
% 
%       - Frequency for amplitude stop (fAend):
%         the end of the examined range of frequencies for amplitude (frequencies of high-frequency oscillation
%         with modulated amplitude). It can not be bigger than sampling frequency/2 divided by 1.15 (which accounts for 
%         filtering transition zone) and with fA_filtband/2 subtracted.
% 
%       - Frequency for amplitude step (fAstep):
%         the step between the beginning and end of examined range of frequencies for amplitude.  
%         
%       - Directory for output (dirOut):
%         directory where the output images and MAT structures will be saved. If there already are the 
%         results of PAC analysis the new results will be saved with added prefix 'vX_', where 'X-1' is the 
%         number of existing sets of results 
% 
%       - Type of data (type):
%         type of EEG data to be processed. It should be 'continuous' when the EEG data is not devided into 
%         epochs or when the existing epochs are supposed to give separate results. It should be 'event_related' 
%         when the epochs are in fact subsequent trials and they are supposed to give one averaged result
% 
%       - Number of phase bins (nbBins):
%         phase of low-frequency oscillation [-pi,pi] will be divided into nbBins bins of equal width. The default value is 18.
%
%       - Wavelet number (number of visible cycles in one burst) (w):
%         number of high-frequency cycles attributable to one transient burst of this frequency. For example, 
%         in the analysis of data related to sequential memory, considering the Lisman model (Lisman, 2005), 
%         it can be set at the level of the number of objects stored in the memory. The default value is set to 5
%
%       - Number of surrogate data repetition (Nboot):
%         number of repetitions while producing the surrogate data. We set the default value to 200
%
%       - Percentile of surrogate data PAC measure distribution (peMI):
%         [1-100] threshold for detecting statistically significant values of eMI. It is expressed as percentile of 
%         distribution of maximal eMI values from each comodulogram for surrogate data.
% 
%       - Percentile of surrogate data Mean Amplitude distribution (pPhaseCom):
%         [1-100] threshold for detecting statistically significant mean amplitudes (averaged for each phase bin). It is expressed as percentile of 
%         distribution of maximal mean amplitude values for surrogate data.
% 
%       - Whether to exclude epochs marked as artifacts (artifacts):
%         [true|false] Epochs marked as artifacts in any of EEG structre fields: EEG.reject.rejjpE, EEG.reject.rejkurtE, 
%         EEG.reject.rejmanualE, EEG.reject.rejthreshE, EEG.reject.rejconstE, EEG.reject.rejfreqE will be excluded from 
%         analysis (only marks for selected channels [chan_fP, chan_fA] are taken into account)
%
%       - Wether to plot with statistical mask (plotWithMask):
%         [1|0] When marked as 0, the generation of surrogate data and statistics will be inactive. When marked as 1, 
%         the surrogate data will be produced, the statistics will be employed and the plots with statistical 
%         mask will be presented
%
%       - Wether to plot without statistical mask (plotWithoutMask):
%         [1|0] When marked as 0, the plots without statistics will not be provided. When marked as 1, the results 
%         without statistics will be presented.
% 
% 
%     Output:
% 
%       All output files are saved in directory specified as dirOut, and {X} below is an identifier of results.
% 
%       - v{X}_eMI_config.mat
%         matlab structure with values of all parameters of PAC analysis described above
% 
%       - v{X}_eMI_results.mat
%         matlab file containing:
%               * fA - vector of frequencies for amplitude
%               * fP - vector of frequencies for phase
%               * fP_bins - array (2 x fP), containing edges of frequency for phase bins
%               * fSpectrum - vector (1 x fAend-fAstart+1), vector of frequencies for amplitude used in averageSpectrum and spectrumOfAverage 
%               * comodulogram - array (fA x fP) with PAC measure for each pair of phase and amplitude frequency
%               * pval - array (fA x fP) with percentiles of surrogate data eMI distibution assigned to each value in comodulogram
%               * comodulogramStat - array (fA x fP), comodulogram with statistical mask on
%               * phaseComodulogram - array (fA x nBins x fP), phase-resolved comodulogram
%               * phaseComodulogramStat - array (fA x nBins x fP), phase-resolved comodulogram with statistical mask on
%               * phaseBins - vector (1 x nBins), phase bins left edges
%               * averageSpectrum - array (fSpectrum x fP), averaged spectrum of signal sections for each fP
%               * spectrumOfAverage - array (fSpectrum x fP), spectrum of averaged signal for each fP 
%               * couplingOrigins - array (fA x fP) indicating origins of coupling: 0-ambiguous or no coupling, 1-reliable coupling
%               * results{fP}Hz - structure with arrays and vectors used to produce the resulting images for each fP:
%                       ** Map - averaged time-frequency map
%                       ** Time - vector of time
%                       ** MeanSignal - averaged signal for amplitude
%                       ** MeanLowFreg - averaged low-frequency oscillation
% 
%       - series of images v{X}_eMI_results{fP}Hz.png (where fP is a frequency from range of frequencies for phase)
%         images containing results of intermediate steps for each fP. The upper plot depicts the averaged 
%         time-frequency map, with black outline of regions that produce a statistically significant coupling.
%         On the right plot there are averaged spectrum (black) and spectrum of averaged signal for amplitude (turquoise). 
%         Shaded area indicates the region of the strongest MI value (grey for ambiguous coupling and blue for reliable coupling). 
%         The blue circle indicates the maximum of spectrum from the range of significant region (black outline) and shaded area.
%         Below you can observe the timecourse of averaged signal for amplitude (turquoise) and averaged low-frequency oscillation (black)
% 
%       - v{X}_eMI_comodulogram.png
%         image of comodulogram without the statistics
% 
%       - v{x}_eMI_comodulogramStat.png
%         image of comodulogram with statistic mask on
%
%       - v{X}_eMI_phaseHistogramStat.png
%         image of polar phase histogram with the statistic mask on. The count values are normalized by the number of all elements in a given region. 
%         The outline color of the phase histogram is consistent with the outline color of the corresponding region in the comodulogramStat.
% 
% 
%     Examples:
%
%     - proper coupling #1: Coupled oscillatory bursts model 
%       more details and description in Test_properPAC.m
%
%     - proper coupling #2: Amplitude modulation model
%       more details and description in Test_properPAC.m
%
%     - proper coupling #3: In vivo recording from rat after ketamine injection 
%       more details and description in Test_properPAC_inVivo.m
%
%     - epiphenomenal coupling #1: Periodic trains of Gaussian functions
%       more details and description in Test_epiphenomenalPAC.m
%
%     - epiphenomenal coupling #2: Non-periodic trains of Gaussian functions
%       more details and description in Test_epiphenomenalPAC.m
%
% -------------------------------------------------------------------------------------------------------------------------

function [LASTCOM] = pop_eMI_calc(EEG , varargin)
    
    LASTCOM = [];
        
    if nargin < 1
        help pop_eMI_calc;
        return;
    end
    
    if isempty(EEG.data)
        error('Cannot process empty dataset');
    end
    
    defaults.chan_fP    = 1;
    defaults.chan_fA    = 1;
    defaults.dirOut     = EEG.filepath;
    defaults.type       = 'continuous';
    defaults.epochs     = 1:EEG.trials; 
    defaults.fPstart    = 4;
    defaults.fPend      = 8;
    defaults.fPstep     = 1;
    defaults.fPband     = 1;
    defaults.fAstart    = 20;
    defaults.fAend      = floor(EEG.srate/2-1);
    defaults.fAstep     = 5;

    defaults.nbBins     = 18;
    defaults.Nboot      = 200;
    defaults.w          = 5;
    defaults.pPhaseCom  = 95;
    defaults.peMI       = 95;
    defaults.artifacts  = true;
    defaults.plotWithMask = 1;
    defaults.plotWithoutMask = 0;
    
    if nargin < 12

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        title_string = 'Detects phase-amplitude coupling by calculating extended Modulation Index -- pop_eMI_calc()';
        geometry = { 1 [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] 1 1 [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] 1 1 [1 1] [1 1] [1 1] [1 1] [1 1]};
        geomvert = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

        uilist = { ...
         { 'style', 'text', 'string', 'Data info: ', 'fontweight', 'bold'}, ...
         { 'style', 'text', 'string', 'Channel number for phase signal: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.chan_fP), 'tag', 'chan_fP'}, ...
         { 'style', 'text', 'string', 'Channel number for amplitude signal: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.chan_fA), 'tag', 'chan_fA'}, ...
         { 'style', 'text', 'string', 'Epoch numbers: ' }, ...
         { 'style', 'edit', 'string', [num2str(min(defaults.epochs)),':',num2str(max(defaults.epochs))], 'tag', 'epochs'}, ...
         { 'style', 'text', 'string', 'Directory for output: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.dirOut), 'tag', 'dirOut'}, ...
         { 'style', 'text', 'string', 'Type of data: ' }, ...
         { 'style', 'popupmenu', 'string', 'continuous|event related' 'tag' 'type' }, ...
         { 'style', 'text', 'string', 'Do you want to analyze epochs marked (if any) as artifacts: ' }, ...
         { 'style', 'popupmenu', 'string', 'no|yes', 'tag', 'artifacts' }, ...
         { 'style', 'text', 'string', 'Which results do you want to plot: ' }, ...
         { 'style', 'checkbox', 'string', {'with statistical mask'}, 'tag', 'plotWithMask', 'Value', 1}, ...
         { },...
         { 'style', 'checkbox', 'string', {'without statistical mask'}, 'tag', 'plotWithoutMask','Value', 0}, ...
         { }, ...
         { 'style', 'text', 'string', 'Analysis parameters: ', 'fontweight', 'bold'}, ...
         { 'style', 'text', 'string', 'Frequency for phase start: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fPstart), 'tag', 'fPstart'}, ...
         { 'style', 'text', 'string', 'Frequency for phase stop: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fPend), 'tag', 'fPend'}, ...
         { 'style', 'text', 'string', 'Frequency for phase step: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fPstep), 'tag', 'fPstep'}, ...
         { 'style', 'text', 'string', 'Frequency for phase band width: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fPband), 'tag', 'fPband'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude start: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAstart), 'tag', 'fAstart'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude stop: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAend), 'tag', 'fAend'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude step: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAstep), 'tag', 'fAstep'}, ...
         { }, ...
         { 'style', 'text', 'string', 'Optional parameters: ', 'fontweight', 'bold'}, ...
         { 'style', 'text', 'string', 'Number of phase bins: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.nbBins), 'tag', 'nbBins'}, ...
         { 'style', 'text', 'string', 'Number of surrogate data repetition: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.Nboot), 'tag', 'Nboot'}, ...
         { 'style', 'text', 'string', 'Wavelet number (number of visible cycles in one burst): ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.w), 'tag', 'w'}, ...
         { 'style', 'text', 'string', 'Percentile of surrogate data MI distribution: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.peMI), 'tag', 'peMI'}, ...
         { 'style', 'text', 'string', 'Percentile of surrogate data Phase Mean Amplitude distribution: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.pPhaseCom), 'tag', 'pPhaseCom'}, ...
         };
     
        [~, ~, err params] = inputgui( 'geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'helpcom', 'pophelp(''pop_eMI_calc'');', 'title' , title_string);
        
        if isempty(params) == 1
            return;
        end
        
        try
            params.chan_fP  = str2num(params.chan_fP);
            params.chan_fA  = str2num(params.chan_fA);
            params.epochs   = str2num(params.epochs);
            params.dirOut   = params.dirOut;
            switch params.type
                case 1
                	params.type = 'continuous';
                case 2
                    params.type = 'event_related';
            end
            switch params.artifacts
                case 1
                	params.artifacts = true;
                case 2
                    params.artifacts = false;
            end
            params.fPstart  = str2num(params.fPstart);
            params.fPend    = str2num(params.fPend);
            params.fPstep   = str2num(params.fPstep);
            params.fPband   = str2num(params.fPband);
            params.fAstart  = str2num(params.fAstart);
            params.fAend    = str2num(params.fAend);
            params.fAstep   = str2num(params.fAstep);

            params.nbBins     = str2num(params.nbBins);
            params.Nboot      = str2num(params.Nboot);
            params.w          = str2num(params.w) ;
            params.peMI       = str2num(params.peMI);
            params.pPhaseCom  = str2num(params.pPhaseCom);

        catch ME1
            throw(ME1);
        end
        
    elseif nargin == 13
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type);
        
        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            
            params.nbBins     = defaults.nbBins;
            params.w          = defaults.w;
            params.Nboot      = defaults.Nboot;
            params.peMI       = defaults.peMI;
            params.pPhaseCom  = defaults.pPhaseCom;
            params.artifacts  = defaults.artifacts;
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
    
    elseif nargin == 14
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            
            params.w          = defaults.nbBins;
            params.Nboot      = defaults.Nboot;
            params.peMI       = defaults.peMI;
            params.pPhaseCom  = defaults.pPhaseCom;
            params.artifacts  = defaults.artifacts;
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
    
    elseif nargin == 15
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};

            params.Nboot      = defaults.Nboot;
            params.peMI       = defaults.peMI;
            params.pPhaseCom  = defaults.pPhaseCom;
            params.artifacts  = defaults.artifacts;
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
        
    elseif nargin == 16
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};
            params.Nboot    = varargin{15};
            
            params.peMI       = defaults.peMI;
            params.pPhaseCom  = defaults.pPhaseCom;
            params.artifacts  = defaults.artifacts;
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
        
    elseif nargin == 17
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};
            params.Nboot    = varargin{15};
            params.peMI     = varargin{16};
            
            params.pPhaseCom  = defaults.pPhaseCom;
            params.artifacts  = defaults.artifacts;
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
        
    elseif nargin == 18
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,pPhaseCom);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};
            params.Nboot    = varargin{15};
            params.peMI     = varargin{16};
            params.pPhaseCom = varargin{17};
            
            params.artifacts  = defaults.artifacts;
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
    
    elseif nargin == 19
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,pPhaseCOm,artifacts);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};
            params.Nboot    = varargin{15};
            params.peMI     = varargin{16};
            params.pPhaseCom = varargin{17};
            params.artifacts = varargin{18};
            
            params.plotWithMask = defaults.plotWithMask;
            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
    
    elseif nargin == 20
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,pPhaseCOm,artifacts,plotWithMask);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};
            params.Nboot    = varargin{15};
            params.peMI     = varargin{16};
            params.pPhaseCom = varargin{17};
            params.artifacts = varargin{18};
            params.plotWithMask = varargin{19};

            params.plotWithoutMask = defaults.plotWithoutMask;
            
        catch ME1
            throw(ME1);
        end
        
    elseif nargin == 21
        % eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbBins,w,Nboot,peMI,pPhaseCOm,artifacts,plotWithMask,plotWithoutMask);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP  = varargin{1};
            params.chan_fA  = varargin{2};
            params.epochs   = varargin{3};
            params.fPstart  = varargin{4};
            params.fPend    = varargin{5};
            params.fPstep   = varargin{6};
            params.fPband   = varargin{7};
            params.fAstart  = varargin{8};
            params.fAend    = varargin{9};
            params.fAstep   = varargin{10};
            params.dirOut   = varargin{11};
            params.type     = varargin{12};
            params.nbBins   = varargin{13};
            params.w        = varargin{14};
            params.Nboot    = varargin{15};
            params.peMI     = varargin{16};
            params.pPhaseCom = varargin{17};
            params.artifacts = varargin{18};
            params.plotWithMask = varargin{19};
            params.plotWithoutMask = defaults.varargin{20};
            
        catch ME1
            throw(ME1);
        end
        
    else
        disp 'Not enough input parameters -- aborting';
        return;
    end
    
    if ((length(size(EEG.data))<3)&&(strcmp(params.type,'event_related')))
        disp('There are no epochs in this dataset. Compulsory change of "type" to continuous')
        params.type = defaults.type;
    end
    
    if (params.fAend>floor(EEG.srate/2))
        disp(['The highest amplitude frequency is the Nyquista frequency. Hence the compulsory change of "Frequency for amplitude stop" to ' num2str(defaults.fAend)])
        params.fAend = defaults.fAend;
    end
    
    if ~(params.plotWithMask || params.plotWithoutMask)
        disp(['You checked not to plot any results. This choice is not optimal, thus the images with statistical mask will be saved.'])
        params.plotWithMask = 1;
    end
    
    if ((params.fPstart-params.fPband/2)<1)
        disp(['The lowest frequency for phase with half of bandwidth subtracted is too close to 0 Hz. Compulsory change of "Frequency for phase start" to ' num2str(round(1+params.fPband/2))])
        params.fPstart = round(1+params.fPband/2);
    end
    
    try
        
        eMI_calc(EEG,params.chan_fP,params.chan_fA,params.epochs,params.fPstart,params.fPend,params.fPstep,params.fPband,params.fAstart,params.fAend,params.fAstep, ...
                 params.dirOut,params.type,params.nbBins,params.w,params.Nboot,params.peMI,params.pPhaseCom,params.artifacts,...
                 params.plotWithMask,params.plotWithoutMask);

             
        tmpstring = '[LASTCOM] = pop_eMI_calc(EEG';
        
        fields = fieldnames(params);
        for ind1 = 1:numel(fields)
            tmpstring = [tmpstring , ' , ' , num2str(params.(fields{ind1}))]; 
        end
        LASTCOM = [tmpstring , ');'];
        
    catch ME1
        throw(ME1);
    end
    
    disp 'Done!'
end
