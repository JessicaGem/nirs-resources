%% May 2021, Jessica Gemignani
%% Generate synthetic fNIRS data

%{

If you use this code in you work, please cite it using one of the following : 
 - Gemignani Jessica, and Judit Gervain (2021) "Comparing different pre-processing routines for infant fNIRS data." Developmental cognitive neuroscience 48: 100943.
 (https://github.com/JessicaGem/mywebsite/blob/gh-pages/GemignaniGervain2021.pdf) %% in this paper the dataset created here is employed to compare pre-processing routines
  
 - Gemignani Jessica (2021) "Classification of fNIRS data with LDA and SVM: a proof-of-concept for application in infant studies" in press in 43rd Annual International Conference of the IEEE Engineering in Medicine and Biology Society
 (pdf available soon) %% in this paper the dataset created here is employed to compare classifiers

- Gemignani Jessica and Gervain Judit (2021) "A practical guide for synthetic fNIRS data generation", in press in 43rd Annual International Conference of the IEEE Engineering in Medicine and Biology Society
(pdf available soon) %% this paper describes each simulation step in detail

%}
%{
Dependencies: 
- Brain AnalyzIR Toolbox (https://github.com/huppertt/nirs-toolbox)
- fixedStimDesign.m function

Inputs: 
- savingFolder
- SD variable (info about montage)
- number of subjects (nsubjects)
- sampling frequency, Hz (Fs)
- tot duration of each timetraces (totSamples)


Outputs: 
- parameters
- data (data_atifacts)

%}

%% To be defined by the user
% savingFolder 
% nsubjects
% SD
% Fs 
% totSamples

%% Setup simulation parameters for three noise levels (low, medium, high)
spike_amplitude=  [2.5 7 18];
shift_amplitude=  [1.5 3 4]; 

subjects_Noise1= nirs.core.Data(); 
subjects_Noise2= nirs.core.Data(); 
subjects_Noise3= nirs.core.Data(); 

parameters_Noise1= cell(1, nsubjects); 
parameters_Noise2= cell(1, nsubjects); 
parameters_Noise3= cell(1, nsubjects); 

for n=1: nsubjects 

    totLength=  totSamples./Fs; % duration 
    t =         (0:1/fs:totLength)'; 
    
    
    %% Step 1) Create experimental design 
    
    % The parameters below replicate the experimental design employed in
    % Gervain et al. 2012. They can be changed as needed
    
    tmin=       300;              % 5 minutes resting state
    stimDur=    18;
    nStimuli=   14;
    stimSpace1= 18+25;            % In this case we want a pause ranging between 25 and 35 s
    stimSpace2= 18+35;
    
    ncond=1;
    
    stim= nirs.testing.fixedStimDesign( t, tmin, stimDur, nStimuli, stimSpace1, stimSpace2, ncond);
    
    %% Step 2) Introduce serial correlations by means of AR model 
    % For more details: Barker et al. 2013
    
    probe= nirs.util.sd2probe(SD);
    
    sigma=              0.33;                % spatial correlation between channels
    ar_order=           30;                % temporal correlation
    
    noise= nirs.testing.simARNoise(probe, t, ar_order, sigma);
    
    %% Step 3) Introduce physiological confound
    % Assume: heart rate ~1.5 Hz, respiration ~0.3 Hz, Meyer waves ~0.1 Hz
    
    cardiacHz=          1.5; 
    respHz=             0.3; 
    meyerHz=            0.1; 
    
    %Specify amplitudes:
    cardiac_amp=        0.25;
    resp_amp=           0.25;
    mayer_amp =         0.25;
    
    data_physio =       nirs.testing.simPhysioNoise_variableHz(noise, cardiac_amp , resp_amp , mayer_amp, ...
                        cardiacHz, respHz, meyerHz);
    
    % TBD: Consolidate ext coefficients and DPF
    
    % ok, now add first HRFs and then motion artifacts (that can ruin HRFs)
    b2 =                0.25;
    b1 =                0.05;
    betaActiveChannels= (b2-b1).*rand(1, 12, 'double'); % created for each subject within the range 0.05 - 0.25;
    
    
    iter=               0;
    
    
    for noiselevel= 1:3
        
        iter=           iter+1;
        spike_amp=      spike_amplitude(noiselevel);
        shift_amp=      shift_amplitude(noiselevel);
        
        
        %% Step 4) Add synthetic HRFs 
        
        % Amplitude of HRF (for HbR it's -HbO/2)

        % HRFs are added only to first 12 channels:
        sd =            unique([noise.probe.link.source noise.probe.link.detector], 'rows');
        channels =      sd(1:round(end/2),:);
        
        % TBD: Consolidate ext coefficients and DPF
        [data, hrf_truth, trueHbO, trueHbR] =     nirs.testing.simData_variableBeta(data_physio, stim, betaActiveChannels, channels );
        
        % Compute initial coefficient of variation, before introducing
        % artifacts
        ini_cv=         nanmean(nanstd(data.data, [], 1)./nanmean(data.data, 1)).*100;
        
        
        %% Step 5) Introduce motion artifacts
        
        spikes_per_min=  2;
        shifts_per_min=  1.5;
        
        % For the typical amplitudes I followed Barker et al 2013
        
        [data_artifacts , truth, spike_inds, shift_inds, baseline_plot_spikes, baseline_plot_shifts] = ...
            nirs.testing.simMotionArtifact(data, spikes_per_minute , shifts_per_minute, spike_amp, shift_amp);
        
        fin_cv= nanmedian(nanstd(data_artifacts.data, [], 1)./nanmean(data_artifacts.data, 1)).*100;
        
        %% Step 6) Save
        
        % Save out useful parameters: 
        
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

        save(fullfile(savingFolder, filename_params), 'parameters', 'ini_cv', 'fin_cv')
        save(fullfile(savingFolder, filename_output), 'data_artifacts')    
            
    end
    
end
