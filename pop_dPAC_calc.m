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


%     dPAC:
%     Calculates the modified Mean Vector Length (direct PAC estimate) - measure of
%     phase-amplitude coupling for given data.
% 
%     Usage:
%       >> pop_dPAC_calc(EEG);          % pop_up window
%
%       >> dPAC_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,...
%                    dirOut,type,Nboot,pdPAC,artifacts,plotWithMask,plotWithoutMask);
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
%       - Filtration bandwidth for phase frequenies (fPfiltband):
%         Each frequency for phase f_P is in fact a bin ranging from f_P-fPfiltband/2 to f_P+fPfiltband/2
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
%       - Filtration bandwidth for amplitude frequency (fAfiltband):
%         Each frequency for amplitude f_A is in fact a bin ranging from f_A-fAfiltband/2 to f_A+fAfiltband/2. 
%         The amplitude-frequency filter width must be at least 2 times bigger than the highest phase-frequency. 
%         It is so, because the spectrum of coupled signal containes side peaks (fA-fP,fA+fP) around coupled high-frequency. 
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
%       - Number of surrogate data repetition (Nboot):
%         number of repetitions while producing the surrogate data. We set the default value to 200
%
%       - Percentile of surrogate data PAC measure distribution (pdPAC):
%         [1-100] threshold for detecting statistically significant values of dPAC. It is expressed as percentile of 
%         distribution of maximal dPAC values from each comodulogram for surrogate data.
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
%       - v{X}_dPAC_config.mat
%         matlab structure with values of all parameters of PAC analysis described above
% 
%       - v{X}_dPAC_results.mat
%         matlab file containing:
%               * fA - vector of frequencies for amplitude
%               * fA_bins - array containing edges of frequency for amplitude bins
%               * fP - vector of frequencies for phase
%               * fP_bins - array containing edges of frequency for phase bins
%               * comodulogram - array (fA x fP) with PAC measure for each pair of phase and amplitude frequency
%               * pval - array (fA x fP) with percentiles of surrogate data dPAC distibution assigned to each value in comodulogram
%               * comodulogramStat - array (fA x fP), comodulogram with statistical mask on

function [LASTCOM] = pop_dPAC_calc(EEG , varargin)

    LASTCOM = [];
        
    if nargin < 1
        help pop_dPAC_calc;
        return;
    end
    
    if isempty(EEG.data)
        error('Cannot process empty dataset');
    end
    
    defaults.chan_fP    = 1;
    defaults.chan_fA    = 1;
    defaults.epochs     = 1:EEG.trials; 
    defaults.dirOut     = EEG.filepath;
    defaults.type       = 'continuous';
    defaults.fPstart    = 4;
    defaults.fPend      = 8;
    defaults.fPstep     = 1;
    defaults.fPfiltband = 1;
    defaults.fAstart    = 20;
    defaults.fAend      = floor(EEG.srate/2/1.15-defaults.fPend);
    defaults.fAstep     = 5;
    defaults.fAfiltband = 2*defaults.fPend;
    defaults.Nboot      = 200;
    defaults.pdPAC      = 95;
    defaults.artifacts  = true;
    defaults.plotWithMask = 1;
    defaults.plotWithoutMask = 0;
    
    if nargin < 18

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        title_string = 'Detects phase-amplitude coupling by calculating direct PAC estimate -- pop_dPAC_calc()';
        geometry = { 1 [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] 1 1 [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1]};
        geomvert = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

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
         { 'style', 'text', 'string', 'Do you want to analyze epochs marked as artifacts: ' }, ...
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
         { 'style', 'text', 'string', 'Frequency for phase filtration band width: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fPfiltband), 'tag', 'fPfiltband'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude start: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAstart), 'tag', 'fAstart'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude stop: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAend), 'tag', 'fAend'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude step: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAstep), 'tag', 'fAstep'}, ...
         { 'style', 'text', 'string', 'Frequency for amplitude filtration band width: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.fAfiltband), 'tag', 'fAfiltband'}, ...
         { 'style', 'text', 'string', 'Number of surrogate data repetition: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.Nboot), 'tag', 'Nboot'}, ...
         { 'style', 'text', 'string', 'Percentile of surrogate data PAC distribution: ' }, ...
         { 'style', 'edit', 'string', num2str(defaults.pdPAC), 'tag', 'pdPAC'}...
         };
     
        [~, ~, err params] = inputgui( 'geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'helpcom', 'pophelp(''pop_dPAC_calc'');', 'title' , title_string);
        
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
            params.fPfiltband = str2num(params.fPfiltband);
            params.fAstart  = str2num(params.fAstart);
            params.fAend    = str2num(params.fAend);
            params.fAstep   = str2num(params.fAstep);
            params.fAfiltband = str2num(params.fAfiltband);
            params.Nboot    = str2num(params.Nboot);
            params.pdPAC    = str2num(params.pdPAC);

        catch ME1
            throw(ME1);
        end
        
    elseif nargin == 19
        % dPAC_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,Nboot,pdPAC,artifacts,plotWithMask,plotWithoutMask);

        if (sum(cellfun('isempty',varargin))>0)
            disp 'There can not be empty parameters -- aborting';
            return;
        end
        
        try
            params.chan_fP    = varargin{1};
            params.chan_fA    = varargin{2};
            params.epochs     = varargin{3};
            params.fPstart    = varargin{4};
            params.fPend      = varargin{5};
            params.fPstep     = varargin{6};
            params.fPfiltband = varargin{7};
            params.fAstart    = varargin{8};
            params.fAend      = varargin{9};
            params.fAstep     = varargin{10};
            params.fAfiltband = varargin{11};
            params.dirOut     = varargin{12};
            params.type       = varargin{13};
            params.Nboot      = varargin{14};
            params.pdPAC      = varargin{15};
            params.artifacts  = varargin{16};
            params.plotWithMask = varargin{17};
            params.plotWithoutMask = varargin{18};
            
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
    
    if (params.fAfiltband < 2*params.fPend)
        disp(['The amplitude-frequency filter width must be at least 2 times bigger than the highest phase-frequency. It is so, because the spectrum of coupled signal containes side peaks (fA-fP,fA+fP) around coupled high-frequency. Hence the compulsory change of "Frequency for amplitude filtration band width:" to ' num2str(2*params.fPend)])
        params.fAfiltband = 2*params.fPend;
    end
    
    if (params.fAstart < params.fPend)
        disp(['The minimal amplitude-frequency is equal to maximal phase-frequency. It is so, because the bandwidth of phase-frequency filter is constant for all amplitude-frequencies and it must be 2 times bigger than the highest phase-frequency. Hence the compulsory change of "Frequency for amplitude start:" to ' num2str(params.fPend)])
        params.fAstart = params.fPend;
    end
    
    if (params.fAend>floor(EEG.srate/2/1.15-params.fAfiltband/2))
        disp(['The highest amplitude frequency is the Nyquista frequency divided by 1.15 (accounts for transision zone), with a half of filtration bandwidth subtracted. Hence the compulsory change of "Frequency for amplitude stop" to ' num2str(floor(EEG.srate/2/1.15-params.fAfiltband/2))])
        params.fAend = floor(EEG.srate/2/1.15-params.fAfiltband/2);
    end
    
    if ~(params.plotWithMask || params.plotWithoutMask)
        disp(['You checked not to plot any results. This choice is not optimal, thus the images with statistical mask will be saved.'])
        params.plotWithMask = 1;
    end
    
    if ((params.fPstart-params.fPfiltband/2)<1)
        disp(['The lowest frequency for phase with half of bandwidth subtracted is too close to 0 Hz. Compulsory change of "Frequency for phase start" to ' num2str(round(1+params.fPfiltband/2))])
        params.fPstart = round(1+params.fPfiltband/2);
    end
    
    try
        
        dPAC_calc(EEG,params.chan_fP,params.chan_fA,params.epochs,params.fPstart,params.fPend,params.fPstep,params.fPfiltband,params.fAstart,params.fAend,...
                       params.fAstep,params.fAfiltband,params.dirOut,params.type,params.Nboot,params.pdPAC,params.artifacts,params.plotWithMask,params.plotWithoutMask);
             
        tmpstring = '[LASTCOM] = pop_dPAC_calc(EEG';
        
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
