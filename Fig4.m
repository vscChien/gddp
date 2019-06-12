% Shih-Cheng Chien, Burkhard Maess, and Thomas R. Knösche. "A generic deviance detection principle for cortical On/Off responses, omission response, and mismatch negativity." BioRxiv (2019): 582437.
function Fig4()
    SOAs=[0 75 125 175 250]; % for Figure 4c
    %SOAs=[0 75:25:250];     % for Figure 4d
    
    
    repetitions=[0 1 2 3 4 5]; % for each SOA
    allR=cell(length(SOAs),length(repetitions));  
    alllatency=zeros(length(SOAs),length(repetitions)); 
    allinput=cell(length(SOAs),length(length(repetitions)));     
    for s=1:length(SOAs) 
       for r=1:length(repetitions) 
          % input: prolonged or periodic stimulus 
          g.soa = SOAs(s);
          g.duration = 50; % msec
          g.onofftime  = [1000 6000-repetitions(r)*max(g.soa,50)]; 
          [input,times]= gen_input(g);
          % inter-column connections
          g.c=[4,2,2,2,1,1,2,2]*0.06;      
          % run simulation
          [~,~,~,~,R]= simulator_21nodes(input,g);
          allR{s,r}=R;
          alllatency(s,r)=find_peak_latency(input,R);
          allinput{s,r}=input;
       end 
    end 
    plot_latency_soa(SOAs,alllatency,g);
    plot_R(allR,allinput,times,SOAs,alllatency,g);
end

%-------------------------------------------------------------------------- 
function [input,times]=gen_input(g)

   dt              = 1e-3;      % 1 msec
   times           = dt:dt:7;   % 7 sec 
   amplitude       = 1.5;       % 1.5 pps
   
   code.n=1;
   code.tlength=length(times);
   code.dt=dt; 
   code.peak=amplitude; 
   code.risetime=10;
   code.offset=g.onofftime(1);
   code.figureon=0;
   if g.soa==0 % prolonged     
       code.soa=20000;
       code.duration=diff(g.onofftime);
       input=gen_input_core(code)'; 
   else        % periodic
       code.soa=g.soa;
       code.duration=g.duration;
       nCycles=floor(diff(g.onofftime)/code.soa);
       input=gen_input_core(code)'; 
       input(g.onofftime(1)+code.soa*nCycles:end)=0; 
   end
     
end
%--------------------------------------------------------------------------  
function plot_latency_soa(SOAs,alllatency,g)
    figure;
    plot(SOAs,alllatency,'o','markeredgecolor',[0 114 189]/255,'markerfacecolor',[0 114 189]/255); hold on;
    plot(SOAs,mean(alllatency,2),'k-o','markeredgecolor',[1 1 1]*0,'markerfacecolor',[1 1 1]*0); 
    plot(SOAs(2:end),SOAs(2:end)-g.duration,'k--'); % due-time 
    xlabel('SOA (ms)');
    ylabel('Peak latency (ms)');
    title('Aligned to offset');
    box on
    ylim([0 600]);xlim([0 max(SOAs)]);
    set(gcf,'position',[0 0 300 300]);
end
%-------------------------------------------------------------------------- 
function plot_R(allR,allinput,times,SOAs,alllatency,g)
    figure;
    for s=1:length(SOAs)
       toff=find(allinput{s,1}~=allinput{s,1}(end),1,'last')+1; 
       subplot(length(SOAs),1,s);hold on
       area(times*1000-toff,(allinput{s,1}>0)*3,'facecolor',[1 1 1]*0.5,'linestyle','none');
       plot(times*1000-toff,allR{s,1},'k-','LineWidth',2)
       if SOAs(s)==0 % due-time 
          plot([0;0],[0;3],'k-','LineWidth',1); 
       else
          plot([1;1]*(SOAs(s)-g.duration),[0;3],'k-','LineWidth',1); 
       end
       scatter(alllatency(s,1),allR{s,1}(toff+alllatency(s,1))+0.5,'v','filled');
       ylim([0 3]); xlim([-800 400]);
       set(gca,'ytick',[],'xtick',-800:200:400,'xminortick','on','ticklength',[0.04 0.025]);
       if SOAs(s)==0
           title('CONT');
       else
           title(['SOA ' num2str(SOAs(s))]);
       end
       set(gca, 'Layer', 'top');
       if s==length(SOAs)
         xlabel('Time (s)');
       end
       box on
    end
    set(gcf,'position',[0 0 300 900]);
end
%-------------------------------------------------------------------------- 
function latency=find_peak_latency(input,R)
   toff=find(input~=input(end),1,'last')+1;               
   [pks,locs] = findpeaks(R(toff+50:toff+500)); 
   [~,maxlocs]=max(pks);
   latency=locs(maxlocs)+50;
end
%--------------------------------------------------------------------------
function data=gen_input_core(code)
   tlength=code.tlength;  
   dt=code.dt;       
   peak=code.peak;      
   soa=code.soa;      
   duration=code.duration;   
   risetime=code.risetime;   
   offset=code.offset;   

   % msec to timepoints
   tlength=tlength/dt/1000;
   soa=soa/dt/1000;
   duration=duration/dt/1000;
   risetime=risetime/dt/1000;
   offset=round(offset/dt/1000);
   % kernel: 50 ms duration (5ms rise and fall)
   u1=interp1([0,risetime,duration-risetime,duration,soa],[0,1,1,0,0]*peak,1:soa); 
   data=repmat(u1,1,ceil(tlength/soa));                        
   data=[zeros(1,offset),data];
   data=data(1:tlength);
   data=data';
end