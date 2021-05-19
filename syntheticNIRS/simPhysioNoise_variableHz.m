function [data, truth] = simPhysioNoise_variableHz( data , cardiac_amp , resp_amp , mayer_amp, cardiacHz, respHz, meyerHz)

% This function is also different from simPhysio in that it does not use
% the DPF and uses different ext coefficients. TBD: consolidate

if nargin<1 || ~exist('data','var') || isempty(data)
   [data, truth]=nirs.testing.simData; 
end
if nargin<2 || isempty(cardiac_amp)
    cardiac_amp = .25; % this gets multiplied by std(hb) to give increase in amplitude
end
if nargin<3 || isempty(resp_amp)
    resp_amp = .25;
end
if nargin<3 || isempty(mayer_amp)
    mayer_amp = .25;
end

if length(data)>1
    for i = 1:length(data)
        data(i) = nirs.testing.simPhysioNoise( data(i) , cardiac_amp , resp_amp , mayer_amp );
    end
    return
end

time = data.time;
[nsamp,nchan] = size(data.data);

% Generate cardiac oscillations
cardiac_freq = cardiacHz + .1*randn; % changed from 1 to 1.5 hz to represent infant's heart rate
cardiac_phase = cumsum( .1*randn(nsamp,1) * 2*pi/data.Fs , 1 );
cardiac_data = sin( 2*pi*cardiac_freq*time + cardiac_phase );

% Generate respiratory oscillations
resp_freq = respHz + .025*randn;
resp_phase = cumsum( .1*randn(nsamp,1) * 2*pi/data.Fs , 1 );
resp_data = sin( 2*pi*resp_freq*time + resp_phase );

% Generate Mayer waves
mayer_freq = meyerHz + .01*randn;
mayer_phase = cumsum( .1*randn(nsamp,1) * 2*pi/data.Fs , 1 );
mayer_data = sin( 2*pi*mayer_freq*time + mayer_phase );

% Add noise to hb concentrations
link = data.probe.link;
Y    = data.data;

if(~(iscellstr(link.type) && any(ismember(link.type,{'hbo','hbr'}))))

    % optical density
    m = mean(Y);
    Y = bsxfun(@plus, -log(Y), log(m));
    
    % sort channels
    channels = nirs.util.uniquerows(link(:,1:2));    
    
    for j = 1:height(channels)
        
        lst = find(link.source==channels.source(j) & link.detector==channels.detector(j));
        
        assert( length(lst) > 1 )

        lambda = link.type(lst);
%         ext = nirs.media.getspectra( lambda );
        e_coef= load('C:\Users\Jessica\Documents\PROJECTS_PARIS\NIRSHitachi_toolbox\NIRSHitachiMehlerPipeline\e_coef.mat'); 
        e_coef= e_coef.e_coef; e_coef= e_coef(:, [1 3 2]); 
        wl1= find(e_coef(:, 1)==lambda(1)); 
        wl2= find(e_coef(:, 1)==lambda(2)); 
        e(1, :)= e_coef(wl1, 2:3); 
        e(2, :)= e_coef(wl2, 2:3); 


        clist = [1 2]; % hbo and hbr; need to fix this

        % extinction coefficients
        E = e; 
        
        % sd distance
        L = data.probe.distances(lst);
        L=max(L,1);  % avoid issues with the short (0) seperation values

        % mbll model
        PPF = 1;
        EL = bsxfun( @times, E, PPF ); %Jessica: removed L and PPF
        iEL = pinv(EL);

        % calculates chromophore concentration (uM)
        hb = (Y(:,lst)*iEL'); %removed 10^-6 factor - these are mM x mm
        
        % add physiological noises
        for k = 1:size(hb,2)
            
            sigma = std(hb(:,k),0,1); % amplitide is increased by cardiac_amp (0.25) * std of the timeseries - 0.0127 mM x mm
            hb(:,k) = hb(:,k) + cardiac_amp * sigma * cardiac_data ...
                              + resp_amp * sigma * resp_data ...
                              + mayer_amp * sigma * mayer_data;
        end
        
        % Convert back to OD
        Y(:,lst) = hb * EL'; %removed 10^-6 factor
        
    end
    
    % Convert OD to intensity
    Y = exp( -bsxfun(@minus, Y, log(m)) );
    
else
    
    % add physiological noise
    for k = 1:size(Y,2)

        sigma = std(Y(:,k),0,1);
        Y(:,k) = Y(:,k) + cardiac_amp * sigma * cardiac_data ...
                        + resp_amp * sigma * resp_data ...
                        + mayer_amp * sigma * mayer_data;
    end

end

data.data = Y;

if(nargout>1)
    if(~exist('truth','var'))
        truth=[];
    end

end