% 20_1 network
%    integrate work of '/media/vchien/V2/_3.code/matlab/20170309_self-sustained_activity (mixed resonance)/demo_RandNetwork_wc_learning.m' 
%    into 'demo4_onoff_periodicinput_wc_20_1.m'
% Periodic input (SOAs=0,75,125,225,300,400 msec)
% layer1:
%    input ----> node1~20
% layer2:
%    (no input)-> node21
% see 6 types of on-off responses at node2
function Fig4()

    %SOAs=[0 75 125 175 250]; % for Figure 4c
    SOAs=[0 75:25:250]; % for Figure 4d
    
    
    repetitions=[0 1 2 3 4 5]; % each SOA
    allR=cell(length(SOAs),length(repetitions));  
    alllatency=zeros(length(SOAs),length(repetitions)); 
    allinput=cell(length(SOAs),length(length(repetitions)));     % for recording
    for s=1:length(SOAs) 
       for r=1:length(repetitions) 
          % Input: prolonged or periodic stimulus 
          g.soa = SOAs(s);
          g.duration = 50; %msec
          g.onofftime  = [1000 6000-repetitions(r)*max(g.soa,50)]; 
          [input,times]=gen_input(g);
          % Inter-column connections
          g.c=[4,2,2,2,1,1,2,2]*0.06;      
          % Run simulation
          [~,~,~,~,R]=simulator_21nodes(input,g);
          allR{s,r}=R;
          alllatency(s,r)=find_peak_latency(input,R);
          allinput{s,r}=input;
       end 
    end 
     plot_latency_soa(SOAs,alllatency,g);
    plot_R(allR,allinput,times,SOAs,alllatency,g);
end

%========================================================================== 
% generate periodic signal
function [input,times]=gen_input(g)

   dt              = 1e-3;     % 1ms
   times           = dt:dt:7; % 16 sec 
   amplitude       = 1.5;      % 1.5 pulses per second
   
   code.type=4; % piecewise
   code.n=1;
   code.tlength=length(times);
   code.dt=dt; 
   code.peak=amplitude; 
   code.risetime=10;
   code.offset=g.onofftime(1);
   code.figureon=0;
   if g.soa==0 % --- square ---      
       code.soa=20000;
       code.duration=diff(g.onofftime);
       input=gen_input_core(code)'; 
   else %-----square periodic----
       code.soa=g.soa;
       code.duration=g.duration;%50;%80;%SOA/2;%ms
       nCycles=floor(diff(g.onofftime)/code.soa);
       input=gen_input_core(code)'; 
       input(g.onofftime(1)+code.soa*nCycles:end)=0; 
   end
     
end
%========================================================================== 
function plot_latency_soa(SOAs,alllatency,g)
    figure;
    plot(SOAs,alllatency,'o','markeredgecolor',[0 114 189]/255,'markerfacecolor',[0 114 189]/255); hold on;
    plot(SOAs,mean(alllatency,2),'k-o','markeredgecolor',[1 1 1]*0,'markerfacecolor',[1 1 1]*0); 
    plot(SOAs(2:end),SOAs(2:end)-g.duration,'k--');% due-time (for reference)
    xlabel('SOA (ms)')
    ylabel('Peak latency (ms)')
    title('Aligned to offset')
    box on
    ylim([0 600]);xlim([0 max(SOAs)]);
    set(gcf,'position',[0 0 300 300])
end
%========================================================================== 
function plot_R(allR,allinput,times,SOAs,alllatency,g)
    figure;
    for s=1:length(SOAs)
       toff=find(allinput{s,1}~=allinput{s,1}(end),1,'last')+1; % time when input turnes off
       subplot(length(SOAs),1,s);hold on
       area(times*1000-toff,(allinput{s,1}>0)*3,'facecolor',[1 1 1]*0.5,'linestyle','none');%green color stimulus
       plot(times*1000-toff,allR{s,1},'k-','LineWidth',2)
       if SOAs(s)==0
          plot([0;0],[0;3],'k-','LineWidth',1); %due-time
       else
          plot([1;1]*(SOAs(s)-g.duration),[0;3],'k-','LineWidth',1); %due-time
       end
       scatter(alllatency(s,1),allR{s,1}(toff+alllatency(s,1))+0.5,'v','filled');
       ylim([0 3])
       xlim([-800 400])
       set(gca,'ytick',[],'xtick',-800:200:400,'xminortick','on','ticklength',[0.04 0.025])
       if SOAs(s)==0
           title('CONT')
       else
           title(['SOA ' num2str(SOAs(s))])
       end
       set(gca, 'Layer', 'top')
       if s==length(SOAs)
         xlabel('Time (s)')
       end
       box on
    end
    set(gcf,'position',[0 0 300 900]);
end
%========================================================================== 
function latency=find_peak_latency(input,R)
   toff=find(input~=input(end),1,'last')+1;  % offset time             
   [pks,locs] = findpeaks(R(toff+50:toff+500)); % find peak in this range [50 500]msec
   [~,maxlocs]=max(pks);
   latency=locs(maxlocs)+50;
end
