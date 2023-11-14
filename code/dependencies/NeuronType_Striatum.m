

function Features=neurontype(DataB)

% input: Data structure of neural recordings for each rat
% output: a matrix called Features: col1: peakvalley, col2: peakwidth,
% col3: neuron type - 0: not classifed, 1: MSN, 2: FSI

Neurons=unique(DataB.Info.NeuronNumb,'rows');
% sampling rate of recordings
sr=24414.0625;

Features=[];

for i=1:size(Neurons,1)%:max(clusters)
timingspikes=(1:size(DataB.Info.Shape(i,:),2))./sr*1000;
Min=min(DataB.Info.Shape(i,:)); %find peaK
idmin=find(DataB.Info.Shape(i,:)==Min);
Max=max(DataB.Info.Shape(i,idmin:end));% find max after peak
idmax=find(DataB.Info.Shape(i,:)==Max);

idhalf=find(DataB.Info.Shape(i,:)<(Min/2));
Peakwidth=(timingspikes(idhalf(end))-timingspikes(idhalf(1))); % Peak width

Peakvalley=(timingspikes(idmax)-timingspikes(idmin)); % Peak valley
Features(i,:)=[Peakvalley Peakwidth 0];   

Type=[];
if Peakvalley>0.560 & Peakvalley<1.5 & Peakwidth>0.15 & Peakwidth<0.45
    Type=1; % Type 1: 
end
    if  Peakvalley>0.1 & Peakvalley<0.5 & Peakwidth>0.05 & Peakwidth<0.2
           Type=2;
    end
    
    if Type==1 | Type==2
        Features(i,3)=[ Type];   

    else
        Type=0;
         Features(i,3)=[Type];   
    end
end
%find(Features(:,3)==0);Features(ans,:)=[];
colors=zeros(size(Features,1),3);
id=find(Features(:,3)==0);colors(id,:)=repmat([0.9290, 0.6940, 0.1250]	,size(id,1),1);
id=find(Features(:,3)==1);colors(id,:)=repmat([1,0,0],size(id,1),1);
id=find(Features(:,3)==2);colors(id,1:3)=repmat([0, 0, 1],size(id,1),1);;

% plot figure S2B
% figure
% scatter(Features(:,1),Features(:,2),120,colors,'filled')
% set(gca,'FontSize',30)
% xlabel('time from peak to valley (ms)')
% ylabel('peak width (ms)')
% axis square
% hold on
% rectangle('Position',[0.56,0.15,0.94,0.3],'EdgeColor','r',...
%     'LineWidth',1)
% hold on
% rectangle('Position',[0.1,0.05,0.4,0.15],'EdgeColor','b',...
%     'LineWidth',1)

%%
end