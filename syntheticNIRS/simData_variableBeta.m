function [data, truth, trueHbO, trueHbR] = simData4Hitachi_variableBeta( noise, stim, beta, channels, basis )
% SIMDATA Simulates NIRS data by adding a task to baseline noise.
%
%     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     October 2019, Jessica Gemignani
%
%     Modifications made:
%     * load extinction coefficients (nirs toolbox has slightly different 
%       values also different order of magnitude (mm-1 x mM-1)
%     * no use of DPF, no use of S-D distances
%     * removed 10^-6 factor
%     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Args:
%     noise  -  a raw data file
%     stim   -  dictionary containing stimulus objects using the stim name as key
%     beta   -  vector of magnitude of response for each stim conditon
%     basis  -  dictionary containing basis objects using stim condition as key
%
% Example:
%     noise = nirs.testing.simARNoise();
%     stim  = nirs.testing.randStimDesign(noise.time, 2, 7, 3);
%     beta  = [3 2 1]';
%
%     sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
%     channels = sd(1:round(end/2),:);
%
%     [data, truth] = simData( noise, stim, beta, channels )


if nargin < 1 || isempty(noise) 
    
    noise = nirs.testing.simARNoise();
end
if strcmp(class(noise),'double')
    noise = nirs.testing.simARNoise(noise);
end

if nargin < 2 || isempty(stim)
    stim = nirs.testing.randStimDesign(noise.time, 2, 7, 1);
end

if nargin < 3 || isempty(beta)
    beta = 7*ones( length(stim.keys), 1 );
elseif(isstr(beta))
    snr = str2num(beta(strfind(beta,'SNR:')+4:end));
    beta=snr*sqrt(var(noise.data(:)));
end

if length(beta) == length(stim.keys)
    % oxy; deoxy
    b = [beta; -beta/2];
else
    b = [beta; -beta/2]; 
end

if nargin < 5 || isempty(basis)
    % default to canonical basis
    basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
end
lstSS=[];

if nargin < 4 || isempty(channels)
    % default to first half of channels
    sd = [noise.probe.link.source noise.probe.link.detector];
    
    % make sure there are no short distance here
    if(ismember('ShortSeperation',noise.probe.link.Properties.VariableNames))
        lstSS=find(noise.probe.link.ShortSeperation);
        sd(lstSS,:)=[];
    else
        lstSS=[];
    end
    sd=unique(sd,'rows');
    channels = sd(1:round(end/2),:);
end

% loop through and add
data = noise.sorted();

link = data.probe.link;
Y    = data.data;
truth = zeros(size(Y,2), 1);

if(~(iscellstr(link.type) && any(ismember(link.type,{'hbo','hbr'}))))
    % optical density
    m = mean(Y);
    Y = bsxfun(@plus, -log(Y), log(m));
    
    
    for i = 1:size(channels, 1)
        lst = find(link.source == channels(i,1) ...
            & link.detector == channels(i,2));
        
        % extincion coefs
        lambda= link.type(lst);
%         e = nirs.media.getspectra(lambda);
%         e = e(:,1:2);
        e_coef= load('C:\Users\Jessica\Documents\PROJECTS_PARIS\NIRSHitachi_toolbox\NIRSHitachiMehlerPipeline\e_coef.mat'); 
        e_coef= e_coef.e_coef; e_coef= e_coef(:, [1 3 2]); 
        
        wl1= find(e_coef(:, 1)==lambda(1)); 
        wl2= find(e_coef(:, 1)==lambda(2)); 
        e(1, :)= e_coef(wl1, 2:3); 
        e(2, :)= e_coef(wl2, 2:3); 
        
        % sd distance
        l = data.probe.distances(lst);
        
        % design mat
        Xhbo = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbo' );
        Xhbr = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbr' );
        
        % entered by Jessica
        Xhbo= sum(Xhbo, 2); 
        Xhbr= sum(Xhbr, 2); 
        
        if size(b, 2)==2
            b=b'; % transpose b for later use (see line 107)
        end
        
        % add to channels according to MBLL
        % changed by Jessica: won't enter s-d distance and DPF (according
        % to later pre-processing)
        
        bCh= b(:, i); 
        
        for j = 1:length(lst)
            % original line
            % Yact = [Xhbo*e(j,1)*l(j) Xhbr*e(j,2)*l(j)] * b * 5/50 * 1e-6;
            % Since this data will be then analyzed not with the nirs
            % toolbox but with the Hitachi pipeline, I remove: DPF, s-d
            % distance and also 10^-6 factor
            Yact = [Xhbo*e(j,1) Xhbr*e(j,2)] * bCh ;

            Y(:,lst(j)) = Y(:,lst(j)) + Yact;
            
            % figure
            %{
            ax= 0:1/10:size(data.data, 1)/10; ax(end)=[];
            figure
            Yact_HbO= Xhbo*e(j,1)*b(1) ;
            Yact_HbO= Xhbr*e(j,2)*b(1) ;
            plot(ax, Yact_HbO)
            xlim([0 50])
            %}
        end
        
        truth(lst) = 1;
    end
    
    Y = exp( -bsxfun(@minus, Y, log(m)) );
else
    for i = 1:size(channels, 1)
        lstHbO = find(link.source == channels(i,1) ...
            & link.detector == channels(i,2) & ismember(link.type,'hbo'));
        lstHbR = find(link.source == channels(i,1) ...
            & link.detector == channels(i,2) & ismember(link.type,'hbr'));
        Xhbo = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbo' );
        Xhbr = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbr' );
        Y(:,lstHbO)=Y(:,lstHbO)+Xhbo*bCh(1);
        Y(:,lstHbR)=Y(:,lstHbR)+Xhbr*bCh(2);
        truth(lstHbO) = 1;
        truth(lstHbR) = 1;
    end
    
    
end

truth(lstSS)=NaN;  % use NaN to mask out the short seperation channels

trueHbO= Xhbo*b(1, :); 
trueHbR= Xhbr*b(2, :); 

data.data = Y;
data.stimulus = stim;
end

