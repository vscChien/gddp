% Shih-Cheng Chien, Burkhard Maess, and Thomas R. Knösche. "A generic deviance detection principle for cortical On/Off responses, omission response, and mismatch negativity." BioRxiv (2019): 582437.
function Fig3()
    % prolonged stimulus (0.5 sec)
    [input,times]=gen_input();
    %---------------Fig. 3e-------------------
    % input connections
    ntrials=100; 
    [ratios1,ratios2]=external_connections_ratios_example1(ntrials);
    % inter-column connections
    g.c=[2 0 2 2 0 0 2 2]*0.1;   % Dec-Off
    % run simulation
    allmu=zeros(ntrials,length(times));
    for i=1:ntrials
        g.ratio1=ratios1(i);
        g.ratio2=ratios2(i);
        g.condition=1; 
        [mu,~,~,~]=simulator_2nodes(input,g);
        allmu(i,:)=mu(2,:); 
    end
    clim=[-0.5 2];
    plot_FRF(allmu,times,ratios1,ratios2,clim);
    
    %---------------Fig. 3f-------------------
    % input connections
    ntrials=100; 
    [ratios1,ratios2]=external_connections_ratios_example2(ntrials);
    % inter-column connections
    g.c=[3 0 2 2 0 0 2 1]*0.1;   % Dec-Off
    % run simulation
    allmu=zeros(ntrials,length(times));
    for i=1:ntrials
        g.ratio1=ratios1(i);
        g.ratio2=ratios2(i);
        g.condition=1; 
        [mu,~,~,~]=simulator_2nodes(input,g);
        allmu(i,:)=mu(2,:); 
    end
    clim=[-0.2 2.5];
    plot_FRF(allmu,times,ratios1,ratios2,clim);
    
    
    %---------------Fig. 3g-------------------
    % input connections
    ntrials=100; 
    [ratios1,ratios2]=external_connections_ratios_example3(ntrials);
    % inter-column connections
    g.c=[4 2 2 2 1 2 2 2 ]*0.1;  % Dec-OnOff
    % run simulation
    allmu=zeros(ntrials,length(times));
    for i=1:ntrials
        g.ratio1=ratios1(i);
        g.ratio2=ratios2(i);
        g.condition=1; 
        [mu,~,~,~]=simulator_2nodes(input,g);
        allmu(i,:)=mu(2,:); 
    end
    clim=[-0.2 2.7];
    plot_FRF(allmu,times,ratios1,ratios2,clim);

end
%--------------------------------------------------------------------------         
function [input,times]=gen_input()

   dt              = 1e-3;        % 1 msec
   times           = dt:dt:5;     % 5 sec 
   amplitude       = 1.5;         % 1.5 pps
   onofftime       = [3000 3500]; % 0.5 sec

   code.n=1;
   code.tlength=length(times);
   code.dt=dt; 
   code.peak=amplitude; 
   code.soa=10000;
   code.duration=diff(onofftime);
   code.risetime=10;
   code.offset=onofftime(1);
   code.figureon=0;
   input=gen_input_core(code)'; 

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
%--------------------------------------------------------------------------  
% ratios on WEX and WIX (for example 1)
function [ratios1,ratios2]=external_connections_ratios_example1(ntrials)
  ratios1=linspace(0.4,1,ntrials);        
  ratios2=gaussmf(linspace(0,1,ntrials),[0.05 0.1])*0.06+...
          gaussmf(linspace(0,1,ntrials),[0.05 0.4])*0.05;       
end
%--------------------------------------------------------------------------  
% ratios on WEX and WIX (for example 2)
function [ratios1,ratios2]=external_connections_ratios_example2(ntrials)
  ratios1=gbellmf(linspace(0,1,ntrials),[0.1 2 0.2])*1;
  ratios2=gbellmf(linspace(0,1,ntrials),[0.05 3 0.4])*.1;   
end
%--------------------------------------------------------------------------  
% ratios on WEX and WIX (for example 3)
function [ratios1,ratios2]=external_connections_ratios_example3(ntrials)
   ratios1=ones(1,ntrials);
   ratios2=gaussmf(linspace(0,1,ntrials),[0.1 0.8])*.45;
end
%--------------------------------------------------------------------------  
function plot_FRF(allmu,times,ratios1,ratios2,clim)
    
    allmu_env=allmu;  % envelope
    for i=1:length(ratios1)
      [pks,locs] = findpeaks(allmu(i,:));
      allmu_env(i,:) = interp1(locs,pks,times*1000);
    end

    figure;
    imagesc(times-3,1:length(ratios1),allmu_env-mean(allmu_env(1,1000:2000)));
    set(gca,'ydir','reverse');
    hold on;
    plot([1;1]*[3 3.5]-3,ylim,'r');
    ylim([0.5 length(ratios1)+0.5]);
    xlim([2.5 4]-3);
    set(gca,'xtick',[0 0.5 1],'xticklabel',[0 0.5 1]);
    xlabel('Time(s)');
    colorbar
    colormap('jet');caxis(clim);
    title('Simulated onset- and offset-FRFs');

    % plot ratios
    pos=get(gca,'position');
    set(gca,'ytick',[],'position',pos+[0.02 0 0 0 ]);
    h=axes('position',[pos(1)-0.08 pos(2) 0.09 pos(4)]);
    hold on
    bar(ratios1,1,'facecolor',[255 102 0]/255,'EdgeColor',[255 102 0]/255);
    bar(ratios2,1,'facecolor','g','EdgeColor','g');
    xlim([0.5 length(ratios1)+0.5]);
    view([-90 -90]);
    box on; set(gca,'layer','top');
    set(gca,'xtick',[],'ytick',[0 1],'yticklabel',{'0' '100%'});
    ylabel('Input ratio');
end

