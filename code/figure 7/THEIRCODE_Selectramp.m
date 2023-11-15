function [Neurons,AllRamps,StereoCrit,MeanFR,preselected_idcs]=...
    THEIRCODE_Selectramp(DataB)
%Function for classify ramping neurons

% input: 1.Data struct of recording session, 2.A NX1 matrix with unique
% NeuronNumber in the session

% extract Neuron Type
Features = THEIRCODE_NeuronType_Striatum(DataB);
Neurons = unique(DataB.Info.NeuronNumb,'rows');

%Pairs of T1-T2
Pairs= [...
    161,334;...
    334,161;...
    334,694;...
    694,334];

% delete the non classifed neurons
idN=find(Features(:,3)==0);
preselected_idcs = Neurons(~ismember(Neurons,idN));
Neurons(idN)=[];

%dt window: size of window for computing PSTH: 75 ms
dt=0.075; %s

% sampling rate of recording session
sr=24414.0625;

%% stereotypy settings
p_fit_idcs = 1 : size(Pairs,1);  % what they say they do in the methods
% p_fit_idcs = 1;                % what's in their code

%% Find ramping neurons aligned to Stimulus 1 onset (window: from -500 to +500 after stim onset)

% preallocation
StereoCrit = struct();
MeanFR = struct();

Ramp_stats=[];Alignment=[];
% loop through all Neurons
for n=1:size(Neurons,1)
    progressreport(n,size(Neurons,1),...
        'flagging ramps and non-ramps (S1 onset)');
    
    % Linear Fit between time and FR
    PSTH=[];Trial_Spikes=[];
    
    %trials from this neuron
    id_N=find(DataB.Info.NeuronNumb==Neurons(n));
    
    for i=1:size(id_N,1)
        Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
            (floor(950):floor(2050))+200),dt*1000,2);
    end
    
    % Calculate PSTH
    PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
    PSTH=PSTH(50:end-50);
    
    % Linear model between time and mean FR
    mdl = fitlm((1:size(PSTH,2))',mean(PSTH,1)');
    % Pearson correlation
    [rho,pval]= corr((1:size(1:size(PSTH,2),2))',PSTH','type','Pearson');
    % polyfit to calculate beta value
    p = polyfit((1:size(PSTH,2))',PSTH',1);
    % concatenate results for all neurons
    Ramp_stats(n,:)=[mdl.Coefficients.pValue(2,1) rho p(1)];
    % Linear Fit between time and FR
    
    Pairs_PSTH=[];Align=[];Trial_Spikes=[];
    for p=1:size(Pairs,1)
        Trial_Spikes=[];
        id_N=find(DataB.Info.NeuronNumb==Neurons(n)&DataB.Info.Duration1==Pairs(p,1));

        if size(id_N,1)>2
            for i=1:size(id_N,1)
                Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
                    (floor(950):floor(2050))+200),dt*1000,2);
            end
            % Calculate PSTH for each T1-T2 pairs
            PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
            Pairs_PSTH(p,:)=PSTH(50:end-50);
        else
            Pairs_PSTH(p,:)=nan(1,size(PSTH(50:end-50),2));
        end
    end
    for p_fit=p_fit_idcs
        Meanvector=mean(Pairs_PSTH,1);
        % corr between each T1-T2 pair PSTH and Mean PSTH
        [R,pval]=corrcoef(Pairs_PSTH(p_fit,:),Meanvector);
        Align(p_fit,:)=[pval(1,2) R(1,2) ];
    end
    
    % find pairs in which correlation is significant and r>0.6
    Stat_Alignement=[Align(:,1)<0.01 Align(:,2)>0.5];
    Alignment(n,:)=[mean(Stat_Alignement(:,1))>=1 mean(Stat_Alignement(:,2))>=1];
    
    %
    StereoCrit.S1_onset(n) = nanmean(Align(:,2));
    MeanFR.S1_onset(n) = nanmean(PSTH);
end

% find significant ramps
id_ramp=find(...
    Ramp_stats(:,1)<0.05&...
    abs(Ramp_stats(:,2))>0.5&...
    abs(Ramp_stats(:,3))>0.004&...
    mean(Alignment,2)==1);

for ii=1:2
    if ii==1
        direction=find(Ramp_stats(id_ramp,3)>0);
        % ramping up neurons
        Stim1onsetup=id_ramp(direction);
    else
        direction=find(Ramp_stats(id_ramp,3)<0);
        % ramping down neurons
        Stim1onsetdown=id_ramp(direction);
    end
end

%% Find ramping neurons aligned to Stimulus 1 offset (window: from -500 to +500 after stim offset)

Ramp_stats=[];Alignment=[];
% loop through all Neurons
for n=1:size(Neurons,1)
    progressreport(n,size(Neurons,1),...
        'flagging ramps and non-ramps (S1 offset)');
    
    % Linear Fit between time and FR
    PSTH=[];Trial_Spikes=[];
    %trials from this neuron
    id_N=find(DataB.Info.NeuronNumb==Neurons(n));
    
    for i=1:size(id_N,1)
        Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
            (floor(950+DataB.Info.Duration1(id_N(i)))...
            :floor(2050+DataB.Info.Duration1(id_N(i))))+200),...
            dt*1000,2);
    end
    % Calculate PSTH
    PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
    PSTH=PSTH(50:end-50);
    
    % Linear model between time and mean FR
    mdl = fitlm((1:size(PSTH,2))',mean(PSTH,1)');
    % Pearson correlation
    [rho,pval]= corr((1:size(1:size(PSTH,2),2))',PSTH','type','Pearson');
    % polyfit to calculate beta value
    p = polyfit((1:size(PSTH,2))',PSTH',1);
    % concatenate results for all neurons
    Ramp_stats(n,:)=[mdl.Coefficients.pValue(2,1) rho p(1)];
    % Linear Fit between time and FR
    
    Pairs_PSTH=[];Align=[];Trial_Spikes=[];
    for p=1:size(Pairs,1)
        Trial_Spikes=[];
        id_N=find(DataB.Info.NeuronNumb==Neurons(n)&DataB.Info.Duration2==Pairs(p,2));

        if size(id_N,1)>2
            for i=1:size(id_N,1)
                Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...,...
                    (floor(950+DataB.Info.Duration1(id_N(i)))...
                    :floor(2050+DataB.Info.Duration1(id_N(i))))+200),...
                    dt*1000,2);
            end
            % Calculate PSTH for each T1-T2 pairs
            PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
            Pairs_PSTH(p,:)=PSTH(50:end-50);
        else
            Pairs_PSTH(p,:)=nan(1,size(PSTH(50:end-50),2));
        end
    end
    for p_fit=p_fit_idcs
        Meanvector=mean(Pairs_PSTH,1);
        % corr between each T1-T2 pair PSTH and Mean PSTH
        [R,pval]=corrcoef(Pairs_PSTH(p_fit,:),Meanvector);
        Align(p_fit,:)=[pval(1,2) R(1,2) ];
    end
    % find pairs in which correlation is significant and r>0.6
    Stat_Alignement=[Align(:,1)<0.01 Align(:,2)>0.5];
    Alignment(n,:)=[mean(Stat_Alignement(:,1))>=1 mean(Stat_Alignement(:,2))>=1];
    
    %
    StereoCrit.S1_offset(n) = nanmean(Align(:,2));
    MeanFR.S1_offset(n) = nanmean(PSTH);
end

% find significant ramps
id_ramp=find(...
    Ramp_stats(:,1)<0.05&...
    abs(Ramp_stats(:,2))>0.5&...
    abs(Ramp_stats(:,3))>0.004&...
    mean(Alignment,2)==1);

for ii=1:2
    if ii==1
        direction=find(Ramp_stats(id_ramp,3)>0);
        % ramping up neurons
        Stim1offsetup=id_ramp(direction);
    else
        direction=find(Ramp_stats(id_ramp,3)<0);
        % ramping down neurons
        Stim1offsetdown=id_ramp(direction);
    end
end

%% Find ramping neurons aligned to Stimulus 2 onset (window: from -500 to +500 after stim onset)

Ramp_stats=[];Alignment=[];
% loop through all Neurons
for n=1:size(Neurons,1)
    progressreport(n,size(Neurons,1),...
        'flagging ramps and non-ramps (S2 onset)');
    
    % Linear Fit between time and FR
    PSTH=[];Trial_Spikes=[];
    %trials from this neuron
    id_N=find(DataB.Info.NeuronNumb==Neurons(n));
    
    for i=1:size(id_N,1)
        Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
            (floor(2950+DataB.Info.Duration1(id_N(i)))...
            :floor(4050+DataB.Info.Duration1(id_N(i))))+200),...
            dt*1000,2);
    end
    % Calculate PSTH for each T1-T2 pairs
    PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
    PSTH=PSTH(50:end-50);
    
    % Linear model between time and mean FR
    mdl = fitlm((1:size(PSTH,2))',mean(PSTH,1)');
    % Pearson correlation
    [rho,pval]= corr((1:size(1:size(PSTH,2),2))',PSTH','type','Pearson');
    % polyfit to calculate beta value
    p = polyfit((1:size(PSTH,2))',PSTH',1);
    % concatenate results for all neurons
    Ramp_stats(n,:)=[mdl.Coefficients.pValue(2,1) rho p(1)];
    % Linear Fit between time and FR
    
    Pairs_PSTH=[];Align=[];Trial_Spikes=[];
    for p=1:size(Pairs,1)
        Trial_Spikes=[];
        id_N=find(DataB.Info.NeuronNumb==Neurons(n)&DataB.Info.Duration2==Pairs(p,2));

        if size(id_N,1)>2
            for i=1:size(id_N,1)
                Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
                    (floor(2950+DataB.Info.Duration1(id_N(i)))...
                    :floor(4050+DataB.Info.Duration1(id_N(i))))+200),...
                    dt*1000,2);
            end
            
            % Calculate PSTH for each T1-T2 pairs
            PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
            Pairs_PSTH(p,:)=PSTH(50:end-50);
        else
            Pairs_PSTH(p,:)=nan(1,size(PSTH(50:end-50),2));
        end
    end
    for p_fit=p_fit_idcs
        Meanvector=mean(Pairs_PSTH,1);
        
        % corr between each T1-T2 pair PSTH and Mean PSTH
        [R,pval]=corrcoef(Pairs_PSTH(p_fit,:),Meanvector);
        Align(p_fit,:)=[pval(1,2) R(1,2) ];
    end
    
    % find pairs in which correlation is significant and r>0.6
    Stat_Alignement=[Align(:,1)<0.01 Align(:,2)>0.5];
    Alignment(n,:)=[mean(Stat_Alignement(:,1))>=1 mean(Stat_Alignement(:,2))>=1];
    
    %
    StereoCrit.S2_onset(n) = nanmean(Align(:,2));
    MeanFR.S2_onset(n) = nanmean(PSTH);
end

% find significant ramps
id_ramp=find(...
    Ramp_stats(:,1)<0.05&...
    abs(Ramp_stats(:,2))>0.5&...
    abs(Ramp_stats(:,3))>0.004&...
    mean(Alignment,2)==1);

for ii=1:2
    if ii==1
        direction=find(Ramp_stats(id_ramp,3)>0);
        % ramping up neurons
        Stim2onsetup=id_ramp(direction);
    else
        direction=find(Ramp_stats(id_ramp,3)<0);
        % ramping down neurons
        Stim2onsetdown=id_ramp(direction);
    end
end

%% Find ramping neurons aligned to Stimulus 2 offset (window: from -500 to +500 after stim offset)

Ramp_stats=[];Alignment=[];
% loop through all Neurons
for n=1:size(Neurons,1)
    progressreport(n,size(Neurons,1),...
        'flagging ramps and non-ramps (S2 offset)');
    
    % Linear Fit between time and FR
    PSTH=[];Trial_Spikes=[];
    %trials from this neuron
    id_N=find(DataB.Info.NeuronNumb==Neurons(n));
    
    for i=1:size(id_N,1)
        Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
            (floor(2950+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i)))...
            :floor(4050+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i))))+200),...
            dt*1000,2);
    end
    % Calculate PSTH for each T1-T2 pairs
    PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
    PSTH=PSTH(50:end-50);
    
    % Linear model between time and mean FR
    mdl = fitlm((1:size(PSTH,2))',mean(PSTH,1)');
    % Pearson correlation
    [rho,pval]= corr((1:size(1:size(PSTH,2),2))',PSTH','type','Pearson');
    % polyfit to calculate beta value
    p = polyfit((1:size(PSTH,2))',PSTH',1);
    % concatenate results for all neurons
    Ramp_stats(n,:)=[mdl.Coefficients.pValue(2,1) rho p(1)];
    % Linear Fit between time and FR
    
    Pairs_PSTH=[];Align=[];Trial_Spikes=[];
    for p=1:size(Pairs,1)
        Trial_Spikes=[];
        id_N=find(DataB.Info.NeuronNumb==Neurons(n)&DataB.Info.Duration2==Pairs(p,2));

        if size(id_N,1)>2
            for i=1:size(id_N,1)
                Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
                    (floor(2950+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i)))...
                    :floor(4050+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i))))+200),...
                    dt*1000,2);
            end
            % Calculate PSTH for each T1-T2 pairs
            PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
            Pairs_PSTH(p,:)=PSTH(50:end-50);
        else
            Pairs_PSTH(p,:)=nan(1,size(PSTH(50:end-50),2));
        end
    end
    for p_fit=p_fit_idcs
        Meanvector=mean(Pairs_PSTH,1);
        % corr between each T1-T2 pair PSTH and Mean PSTH
        [R,pval]=corrcoef(Pairs_PSTH(p_fit,:),Meanvector);
        Align(p_fit,:)=[pval(1,2) R(1,2) ];
    end
    % find pairs in which correlation is significant and r>0.6
    Stat_Alignement=[Align(:,1)<0.01 Align(:,2)>0.5];
    Alignment(n,:)=[mean(Stat_Alignement(:,1))>=1 mean(Stat_Alignement(:,2))>=1];
    
    %
    StereoCrit.S2_offset(n) = nanmean(Align(:,2));
    MeanFR.S2_offset(n) = nanmean(PSTH);
end

% find significant ramps
id_ramp=find(...
    Ramp_stats(:,1)<0.05&...
    abs(Ramp_stats(:,2))>0.5&...
    abs(Ramp_stats(:,3))>0.004&...
    mean(Alignment,2)==1);

for ii=1:2
    
    if ii==1
        direction=find(Ramp_stats(id_ramp,3)>0);
        % ramping up neurons
        Stim2offsetup=id_ramp(direction);
    else
        direction=find(Ramp_stats(id_ramp,3)<0);
        % ramping down neurons
        Stim2offsetdown=id_ramp(direction);
    end
end

%% Find ramping neurons aligned to Go Cue (window: from -500 to +500 after go cue)

Ramp_stats=[];Alignment=[];
% loop through all Neurons
for n=1:size(Neurons,1)
    progressreport(n,size(Neurons,1),...
        'flagging ramps and non-ramps (go cue)');
    
    % Linear Fit between time and FR
    PSTH=[];Trial_Spikes=[];
    %trials from this neuron
    id_N=find(DataB.Info.NeuronNumb==Neurons(n));
    
    for i=1:size(id_N,1)
        Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
            (floor(3450+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i)))...
            :floor(4550+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i))))+200),...
            dt*1000,2);
    end
    % Calculate PSTH for each T1-T2 pairs
    PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
    PSTH=PSTH(50:end-50);
    
    % Linear model between time and mean FR
    mdl = fitlm((1:size(PSTH,2))',mean(PSTH,1)');
    % Pearson correlation
    [rho,pval]= corr((1:size(1:size(PSTH,2),2))',PSTH','type','Pearson');
    % polyfit to calculate beta value
    p = polyfit((1:size(PSTH,2))',PSTH',1);
    % concatenate results for all neurons
    Ramp_stats(n,:)=[mdl.Coefficients.pValue(2,1) rho p(1)];
    % Linear Fit between time and FR
    
    Pairs_PSTH=[];Align=[];Trial_Spikes=[];
    for p=1:size(Pairs,1)
        Trial_Spikes=[];
        id_N=find(DataB.Info.NeuronNumb==Neurons(n)&DataB.Info.Duration2==Pairs(p,2));
        if size(id_N,1)>2
            for i=1:size(id_N,1)
                Trial_Spikes(i,:)=movsum(DataB.Info.FR(id_N(i),...
                    (floor(3450+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i)))...
                    :floor(4550+DataB.Info.Duration1(id_N(i))+DataB.Info.Duration2(id_N(i))))+200),...
                    dt*1000,2);
            end
            % Calculate PSTH for each T1-T2 pairs
            PSTH=nansum(Trial_Spikes,1)./(size(id_N,1)*dt);
            Pairs_PSTH(p,:)=PSTH(50:end-50);
        else
            Pairs_PSTH(p,:)=nan(1,size(PSTH(50:end-50),2));
        end
    end
    for p_fit=p_fit_idcs
        Meanvector=mean(Pairs_PSTH,1);
        % corr between each T1-T2 pair PSTH and Mean PSTH
        [R,pval]=corrcoef(Pairs_PSTH(p_fit,:),Meanvector);
        Align(p_fit,:)=[pval(1,2) R(1,2) ];
    end
    % find pairs in which correlation is significant and r>0.6
    Stat_Alignement=[Align(:,1)<0.01 Align(:,2)>0.5];
    Alignment(n,:)=[mean(Stat_Alignement(:,1))>=1 mean(Stat_Alignement(:,2))>=1];
    
    %
    StereoCrit.Go_cue(n) = nanmean(Align(:,2));
    MeanFR.Go_cue(n) = nanmean(PSTH);
end

% find significant ramps
id_ramp=find(...
    Ramp_stats(:,1)<0.05&...
    abs(Ramp_stats(:,2))>0.5&...
    abs(Ramp_stats(:,3))>0.004&...
    mean(Alignment,2)==1);

for ii=1:2
    if ii==1
        direction=find(Ramp_stats(id_ramp,3)>0);
        % ramping up neurons
        gocueup=id_ramp(direction);
    else
        direction=find(Ramp_stats(id_ramp,3)<0);
        % ramping down neurons
        gocuedown=id_ramp(direction);
    end
end

%%
%Create an output matrix called AllRamps with N rows for each Neuron, and
%10 colums. Each colum indicate a type or ramping neuron:
%1) Aligned to Stim1 onset - Ramping UP,
%2) Aligned to Stim1 onset - Ramping DOWN
%3) Aligned to Stim1 offset - Ramping UP,
%4) Aligned to Stim1 offset - Ramping DOWN
%5) Aligned to Stim2 onset - Ramping UP,
%6) Aligned to Stim2 onset - Ramping DOWN
%7) Aligned to Stim2 offset - Ramping UP,
%8) Aligned to Stim2 offset - Ramping DOWN
%9) Aligned to GoCue  - Ramping UP,
%10) Aligned to GoCue  - Ramping DOWN

%entries are either 1 or 0, indicated if the Neuron at each Row is classifed as
%any of the 10 type of ramps
AllRamps=zeros(size(Neurons,1),10);
AllRamps(Stim1onsetup,1)=1;
AllRamps(Stim1onsetdown,2)=1;
AllRamps(Stim1offsetup,3)=1;
AllRamps(Stim1offsetdown,4)=1;
AllRamps(Stim2onsetup,5)=1;
AllRamps(Stim2onsetdown,6)=1;
AllRamps(Stim2offsetup,7)=1;
AllRamps(Stim2offsetdown,8)=1;
AllRamps(gocueup,9)=1;
AllRamps(gocuedown,10)=1;
end