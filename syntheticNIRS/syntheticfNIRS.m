%% Jan 2021, Jessica Gemignani
%% Generate synthetic dataset

% Description: In order to generate a synthetic dataset based on a real one, a "bb"
% structure from previous analysis on real data is needed; it will be used
% to extract name, stimulus design, sampling frequency
% Additionally, the SD (probe arrangement) is necessary, to replicate it

% Output: this code takes the csv file of a subject, copies it to another location
% and overwrites the real data with the synthetic data (maintaining
% stimulus design and everything else intact)

% *** Left to be done ***: also save out characteristics of 
% - HRF amplitude
% - AR model of the baseline noise
% - physiological confound characteristics
% - motion artifacts' locations and amplitudes 
%
% all of these might be useful when evaluating performance of different methods

nirstoolbox= 'C:\Users\Jessica\Documents\Matlab_Toolboxes\nirs-toolbox-clone2'; 
addpath(genpath(nirstoolbox))

load('bb_rr.mat')
nsubjects= numel(bb); 
load('SD_new.mat')
% 10 subjects
bb(11:end)=[]; 

%% Select subject

% Synthetic data will be generated based on this subject experimental
% design. 
% The new data will be saved in a csv file with the same name as
% the original file + "_simulated". Eventually, many iterations can be
% performed on the same subject - to produce a robust statistics - each
% time tuning the HRF, physiological confound and motion artifacts
% parameters.

%% Setup simulation parameters

HRF_amp =         [0.1 0.2 0.3];
spike_amplitude=  [2.5 7 18];
shift_amplitude=  [1.5 3 4]; 

%%
savingFolder= 'C:\Users\Jessica\Documents\PROJECTS_PADOVA\2021_LDA_SVM_EMBC\SyntheticData';

subjects_Noise1= nirs.core.Data(); 
subjects_Noise2= nirs.core.Data(); 
subjects_Noise3= nirs.core.Data; 

parameters_Noise1= cell(1, 10); 
parameters_Noise2= cell(1, 10); 
parameters_Noise3= cell(1, 10); 

for n=1: 10% nsubjects % subject n1
    
    totSamples= size(bb(n).data, 1);
    fs=         1./bb(n).measure.samplingperiod;
    totLength=  totSamples./fs; %length of measurement in seconds
    t = (0:1/fs:totLength)'; % time vector
    
    
    %% Create experimental design of that subject
    %     stim= extractExperimentalDesign(bb(n));
    % then exclude N and make only R trials
    
    % Visualize stimulus
    %{
                t=1:totLength;
                vec1= stim.values{1}.getStimVector(t); vec2= stim.values{2}.getStimVector(t);
                figure, plot(ax, vec1, 'k', 'LineWidth', 1.5),
                %hold on, plot(t, vec2, 'k', 'LineWidth', 1.5)
                ylim([0 1]), %legend('ABC', 'AAB')
                axis off
                set(gca, 'LineWidth', 2)
    %}
    
    tmin= 300; % 5 minutes resting state
    stimDur= 18;
    nStimuli= 14;
    stimSpace1= 18+25;
    stimSpace2= 18+35;
    ncond=1;
    stim= nirs.testing.fixedStimDesign( t, tmin, stimDur, nStimuli, stimSpace1, stimSpace2, ncond);
    
    %% Step 1) AR noise + physiological confound
    probe= nirs.util.sd2probe(SD);
    
    sigma= 0.33; % spatial correlation between channels
    ar_order=30; % temporal correlation
    
    noise= nirs.testing.simARNoise(probe, t, ar_order, sigma);
    
    %% Step 1) Introduce physiological confound
    % Assume: heart rate ~1.5 Hz, respiration ~0.3 Hz, Meyer waves ~0.1 Hz
    
    %Specify amplitudes:
    cardiac_amp= 0.25;
    resp_amp= 0.25;
    mayer_amp = 0.25;
    
    data_physio =   nirs.testing.simPhysioNoise4Hitachi(noise, cardiac_amp , resp_amp , mayer_amp);
    
    % ok, now add first HRFs and then motion artifacts (that can ruin HRFs)
    b2 = 0.25;
    b1 = 0.05;
    betaActiveChannels= (b2-b1).*rand(1, 12, 'double'); % created for each subject within the range 0.05 - 0.25;
    
    
    iter=0;
    
    
    for noiselevel= 1:3
        
        iter=iter+1;
        spike_amp= spike_amplitude(noiselevel);
        shift_amp= shift_amplitude(noiselevel);
        
        
        %% Step 2) Add synthetic HRFs based on extracted experimental design
        
        % Amplitude of HRF (for HbR it's -HbO/2)
        % beta=   0.2;
        % HRFs are added only to first 12 channels:
        sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
        channels = sd(1:round(end/2),:);
        
        [data, hrf_truth, trueHbO, trueHbR] =     nirs.testing.simData4Hitachi_VariableBeta(data_physio, stim, betaActiveChannels, channels );
        
        ini_cv= nanmean(nanstd(data.data, [], 1)./nanmean(data.data, 1)).*100;
        %{
%                 figure
%                 plot(ax, data.data(:, 1), 'Color', grey)
%                 ylim([0.9 1.1])
%                 xlim([0 150])
%                 xlabel('seconds'), ylabel('V')
%                 set(gca, 'FontSize', 18)
% %                 saveas(gcf, fullfile('.\Figures', 'dataHRF_rawData.fig'), 'fig')
%
%                 hb_HRF= job.run(data);
%
%                 figure
%                 plot(ax, hb_HRF.data(:, 1), 'Color', grey)
% %                 ylim([0.9 1.1])
%                 xlim([0 150])
%                 xlabel('seconds'), ylabel('mM x mm')
%                 set(gca, 'FontSize', 18)
% %                 saveas(gcf, fullfile('.\Figures', 'dataHRF_hbo.fig'), 'fig')
%
%                 figure
%                 plot(ax, hb_HRF.data(:, 1), 'Color', 'r')
%                 hold on
%                 plot(ax, hb_HRF.data(:, 2), 'Color', 'b')
%                 xlabel('seconds'), ylabel('mM x mm')
%                 set(gca, 'FontSize', 22)
% %                 saveas(gcf, fullfile('.\Figures', 'dataHRF_hbo_hbr.fig'), 'fig')
        %}
        
        %% Step 3) Introduce motion artifacts
        
        spikes_per_minute=  2;
        shifts_per_minute=  1.5;
        
        % For the typical amplitudes I followed Barker et al 2013
        
        [data_artifacts , truth, spike_inds, shift_inds, baseline_plot_spikes, baseline_plot_shifts] = ...
            nirs.testing.simMotionArtifact(data, spikes_per_minute , shifts_per_minute, spike_amp, shift_amp);
        
        fin_cv= nanmedian(nanstd(data_artifacts.data, [], 1)./nanmean(data_artifacts.data, 1)).*100;
        
        %{
%                 figure, plot(data_artifacts.data)
%                 job= nirs.modules.OpticalDensity();
%                 job= nirs.modules.BeerLambertLaw4Hitachi(job);
%                 hb_art= job.run(data_artifacts);
%                 figure, plot(hb_art.data)
%
%
%                 figure,
%                 plot(ax, baseline_plot_shifts, 'Color', rgb('DarkGreen'), 'LineWidth', 3)
%                 hold on
%                 plot(ax, baseline_plot_spikes,'Color', rgb('Purple'), 'LineWidth', 3)
%
%                 figure
%                 plot(ax, data_artifacts.data(:, 1), 'Color', grey)
% %                 ylim([0.9 1.1])
% %                 xlim([0 500])
%                 xlabel('seconds'), ylabel('V')
%                 set(gca, 'FontSize', 30)
% %                 saveas(gcf, fullfile('.\Figures', 'dataArtifacts_rawData.fig'), 'fig')
%
%                 job= nirs.modules.OpticalDensity();
%                 job= nirs.modules.BeerLambertLaw4Hitachi(job);
%                 hb_unf= job.run(data_artifacts);
%
%                 figure, plot(ax, hb_unf.data(:, 1), 'r', 'LineWidth', 1.2)
%                 hold on, plot(ax, hb_unf.data(:, 2), 'b', 'LineWidth',1.2)
%                 ylim([-2.5 2.5]),
%                 xlim([0 1500])
%                 xlabel('seconds'), ylabel('mM x mm')
%                 set(gca, 'FontSize', 40)
%                 saveas(gcf, fullfile('.\Figures', 'dataArtifacts_hb.fig'), 'fig')

%% Quantify the magnitude of signal change introduced by the artifacts
%                 yy= data_artifacts.data(:, 1);
%                 disp('min:')
%                 min(yy)
%                 disp('max:')
%                 max(yy)
%                 disp('mean:')
%                 mean(yy)
                
%                 figure,
%                 plot(yy)
%                 hold on
%                 plot([spike_inds' spike_inds'], get(gca, 'ylim'), 'r')
%                 hold on
%                 plot([(spike_inds-100)' (spike_inds-100)'], get(gca, 'ylim'), 'k')
%                                hold on
%                 plot([(spike_inds+100)' (spike_inds+100)'], get(gca, 'ylim'), 'k')
%
%                 for inds=1:numel(spike_inds)
%                     if spike_inds>100
%                     peakMax(inds) = max(yy(spike_inds(inds)-100:spike_inds(inds)+100));
%                     peakMin(inds) = min(yy(spike_inds(inds)-100:spike_inds(inds)+100));
%                     else continue
%                     end
%                 end
%
%                 spike_changes= peakMax-peakMin;
%                 disp('changes for spikes:')
%                 min(spike_changes)
%                 max(spike_changes)
%                 max(peakMax)
%
%                 % shift changes
%
%                 % quantify magnitude of artifacts
%                 yy= data_artifacts.data(:, 1);
%                 figure,
%                 plot(yy)
%                 hold on
%                 plot([shift_inds' shift_inds'], get(gca, 'ylim'), 'r')
%                 hold on
%                 plot([(shift_inds-25)' (shift_inds-25)'], get(gca, 'ylim'), 'k')
%                                hold on
%                 plot([(shift_inds+25)' (shift_inds+25)'], get(gca, 'ylim'), 'k')
%
%                 for inds=1:numel(shift_inds)
%                     peakMax(inds) = max(yy(shift_inds(inds):shift_inds(inds)+25));
%                     peakMin(inds) = min(yy(shift_inds(inds)-25:shift_inds(inds)));
%                 end
%
%                 shift_changes= peakMax-peakMin;
%                 disp('changes for shifts:')
%                 min(shift_changes)
%                 max(shift_changes)
%                 max(peakMax)

                
%                 figure
%                 plot(ax, data_artifacts.data(:, 1), 'Color', grey)
%                 ylim([0 2.5])
%                 xlim([0 1400])
%                 xlabel('seconds'), ylabel('V')
%                 set(gca, 'FontSize', 18)
                
        %{
                %job= nirs.modules.OpticalDensity(); job= nirs.modules.BeerLambertLaw4Hitachi(job);
                % hb= job.run(data);
                % hb_artifacts= job.run(data_artifacts);
                %
                % figure,
                % plot(hb.data(:, 1))
                % hold on
                % plot(hb_artifacts.data(:, 1), 'r')
                % title('spike amp= 3, shift amp= 5')
                % ylim([-4 10])
        %}

                % See generated data:
        %{
                figure,
                plot(noise.data(:, 1), 'k')
                hold on
                plot(data_physio.data(:, 1), 'b')
                hold on
                plot(data.data(:,1), 'r')
                hold on
                plot(data_artifacts.data(:, 1), 'm')
                ylim([0 10])
                legend('baseline', 'baseline + physio', 'bl + physio + HRF', '+artifacts')
                hold on
                plot([shift_inds' shift_inds'], get(gca, 'ylim'), 'k')
        %}

%                 %% Now data is ready (and stored in "data_artifacts")
%                 % Write data into csv file
%                 filename_input= [bb(n).id '.csv'];
%                 filename_output= [bb(n).id '_' iteration_name '.csv'];
%
%                 inputFolder= 'C:\Users\Jessica\Documents\PROJECTS_PARIS\InfantsPreprocessingPipelines_code\RealDataBoysGirls\PipelineA\RawData';
%                 outputFolder= 'C:\Users\Jessica\Documents\PROJECTS_PARIS\InfantsPreprocessingPipelines_code\SimulatedData\GeneratedData';
%
%                 % Pick CSV file to modify and copy it to destination folder first before
%                 % overwriting data
%                 copyfile (fullfile(inputFolder, filename_input), fullfile(outputFolder, filename_output))
%
%                 % data to be written
%                 A= data_artifacts.data;
%                 A(end, :)=[]; %for some reason the simulation adds one extra value (one more time sample)
%
%                 % data in the csv file starts at row 42
%                 startRow= 42;
%                 endRow= size(A, 1) + startRow  -1;
%                 range1= ['B' num2str(startRow)];
%                 range2= ['AW' num2str(endRow)];
%                 range= [range1 ':' range2];

                % write data
%                 xlswrite(fullfile(outputFolder, filename_output),A, range)
        %}
        
        % I also want to save out parameters related to this
        % particular subject and iteration
        parameters.iteration=   iter;
        parameters.subject=     n;
        parameters.beta=        betaActiveChannels;
        parameters.active_ch=   hrf_truth;
        parameters.trueHbO=     trueHbO;
        parameters.trueHbR=     trueHbR;
        parameters.spike_inds=  spike_inds;
        parameters.shift_inds=  shift_inds;
        parameters.spike_amp=   spike_amp;
        parameters.shift_amp=   shift_amp;
        parameters.spike_min=   spikes_per_minute;
        parameters.shift_min=   shifts_per_minute;
        parameters.cardiac_amp= cardiac_amp;
        parameters.resp_amp=    resp_amp;
        parameters.mayer_amp=   mayer_amp;
        parameters.ar_order=    ar_order;
        
        if noiselevel==1
            subjects_Noise1(n)=         data_artifacts; 
            parameters_Noise1{1, n}=    parameters; 
        elseif noiselevel==2
            subjects_Noise2(n)=         data_artifacts; 
            parameters_Noise2{1, n}=    parameters; 
        elseif noiselevel==3
            subjects_Noise3(n)=         data_artifacts; 
            parameters_Noise3{1, n}=    parameters; 
        end
        
%         filename_output= ['Subject' num2str(n) '_VariableBeta_NoiseLevel' num2str(noiselevel) '.mat'];
%         filename_params= [filename_output(1:end-4) '_Parameters.mat'];
%         save(fullfile(savingFolder, filename_params), 'parameters')
%         save(fullfile(savingFolder, filename_output), 'data_artifacts')
        
        
        %{
                vec1= stim.values{1}.getStimVector(t); vec2= stim.values{2}.getStimVector(t);

                %% Processing: Optical Density and Beer Lambert
                job= nirs.modules.OpticalDensity();
                job= nirs.modules.BeerLambertLaw4Hitachi(job);
                job.PPF=[1 1];

                hb= job.run(data);
                hb_physio= job.run(data_physio);
                hb_artif= job.run(data_artifacts);
                hb_noise= job.run(noise);

                %% Bandpass filter
                job= eeg.modules.BandPassFilter();
                job.do_downsample=0;
                job.lowpass=0.5;
                job.highpass=0.01;

                filtNoise=      job.run(noise);
                filtHb=         job.run(hb);
                filtHb_physio=  job.run(hb_physio);
                filtHb_artif=   job.run(hb_artif);

                %% Plots
                figure,
                subplot(2,1,1)
                plot(filtNoise.data(:, 1), 'k', 'LineWidth', 1.5)
                hold on
                plot(filtHb.data(:, 1), 'r', 'LineWidth', 1.5)
                hold on
                plot(filtHb_physio.data(:, 1), 'm', 'LineWidth', 1.5)
                hold on
                plot(filtHb_artif.data(:, 1), 'b', 'LineWidth', 1.5)
                hold on
                grey=[130 130 130]./255;
                plot(vec1+vec2, 'Color', grey', 'LineWidth',2)
                legend('baseline (bl)', 'bl+HRF', 'bl+HRF+physiological', 'bl+HRF+physio+artifacts', 'stim design')
                title('SIMULATION (Band-pass filtered HbO)')
                ylim([-5 5])

                subplot(2,1,2)
                plot(bb(n).hb(:, 1,1), 'b', 'LineWidth', 1.5)
                hold on
                grey=[130 130 130]./255;
                plot(vec1+vec2, 'Color', grey', 'LineWidth',2)
                legend('HbO', 'stim design')
                title('REAL MEASURE (Band-pass filtered HbO)')
                ylim([-5 5])

                figure,
                plot(hb_noise.data(:, 1), 'k', 'LineWidth', 1.5)
                hold on
                plot(hb.data(:, 1), 'r', 'LineWidth', 1.5)
                hold on
                plot(hb_physio.data(:, 1), 'm', 'LineWidth', 1.5)
                hold on
                plot(hb_artif.data(:, 1), 'b', 'LineWidth', 1.5)
                hold on
                grey=[130 130 130]./255;
                plot(vec1+vec2, 'Color', grey', 'LineWidth',2)
                legend('baseline (bl)', 'bl+HRF', 'bl+HRF+physiological', 'bl+HRF+physio+artifacts', 'stim design')
                title('SIMULATION (unfiltered HbO)')
                ylim([-5 5])



        %}
        
    end
    
    % before passing onto the next iteration write out anything useful about
    % this one
    
end