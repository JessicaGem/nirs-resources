function stim = fixedStimDesign( t, tmin, stimDur, nStimuli, stimSpaceLB, stimSpaceUB, ncond, condnames)

%% fixedStimDesign - Returns a stim design with predefined timing of stimuli and pauses

%{

Based on nirs.testing.randStimDesign() (BrainAnalyzIR Toolbox), this function allows to predefine 
the number of stimuli to be simulated and their timings, as well as the
timings of the interstimuli intervals. 

-------------------------
Jessica Gemignani

% Args:
%     t             - time vector (s)
%     tmin          - initial segment of resting state, duration (s)
%     stimDur       - stim duration (s)
%     nStimuli      - tot number of trials in the whole experiment (i.e: num trials per condition = nStimuli/ncond)
%     stimSpaceLB   - minimum space between stim onsets (s)
%     stimSpaceUB   - maximum space between stim onsets (s)
%     ncond         - number of conditions
%     condnames     - cell-array with names of conditions, e.g. {'Play', 'Rest'}


%}

    if nargin==8 & length(condnames) ~= ncond
        disp('Attention: names and number of conditions do not match')
        return
    end

    if nargin < 8 %missing condnames
        condnames= cell(1, ncond); 
        disp('Conditions will be automatically named')
        for ii=1:ncond
            condnames{ii}= ['condition_' num2str(ii)]; 
            
        end
    end
        
        
    if nargin < 7 % missing number of conditions
        ncond = 1;
    end
    
    % min max times
    %tmin = min( t ) + 1*stimDur; % change here the tmin
    tmax = max( t ) - 2*stimDur;

    % number of possible stims onsets
    nrnd = round( 2*(tmax-tmin)/stimSpaceUB );
    
    % random times between tasks
    % dt = stimSpace/2 + exprnd(stimSpace/2, [nrnd 1]); 
    dt= randi([stimSpaceLB stimSpaceUB], nrnd, 1); 

    % onsets
    onset = tmin + cumsum([0; dt]);
    onset = onset( onset < tmax );
    % now keep only nStimuli
    onset= onset(randperm(numel(onset), nStimuli)); 
    onset= sort(onset); 
    

    % durations
    dur = stimDur * ones(size(onset));

    % amplitude
    amp = ones(size(dur));

    % output
    stim = Dictionary();
    r = randi(ncond);
    
    for i = 1:ncond
        lst = mod(i+r,ncond)+1:ncond:size(dur,1);     
        s = nirs.design.StimulusEvents();
        s.name   = condnames{i};
        s.amp    = amp(lst);
        s.dur    = dur(lst);
        s.onset  = onset(lst);
        stim(s.name) = s;
    end
    
end

