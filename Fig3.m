% Comparison Between Offset and Onset Responses of Primary Auditory Cortex on–off Neurons in Awake Cats, 2007
function Fig3()

    % Input: prolonged stimulus
    [input,times]=gen_input();
    %---------------Fig. 3e-------------------
    % Input connections
    ntrials=100; 
    [ratios1,ratios2]=external_connections_ratios_example1(ntrials);
    % Inter-column connections
    g.c=[2 0 2 2 0 0 2 2]*0.1;   % inhibited-Off
    % Run simulation
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
    % Input connections
    ntrials=100; 
    [ratios1,ratios2]=external_connections_ratios_example2(ntrials);
    % Inter-column connections
    g.c=[3 0 2 2 0 0 2 1]*0.1;   % inhibited-Off
    % Run simulation
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
    % Input connections
    ntrials=100; 
    [ratios1,ratios2]=external_connections_ratios_example3(ntrials);
    % Inter-column connections
    g.c=[4 2 2 2 1 2 2 2 ]*0.1;  % inhibited-On&Off
    % Run simulation
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
%==========================================================================
% ratios on WEX and WIX
function [ratios1,ratios2]=external_connections_ratios_example1(ntrials)
  ratios1=linspace(0.4,1,ntrials);        
  ratios2=gaussmf(linspace(0,1,ntrials),[0.05 0.1])*0.06+...
          gaussmf(linspace(0,1,ntrials),[0.05 0.4])*0.05;       
end
%==========================================================================
% ratios on WEX and WIX
function [ratios1,ratios2]=external_connections_ratios_example2(ntrials)
  ratios1=gbellmf(linspace(0,1,ntrials),[0.1 2 0.2])*1;
  ratios2=gbellmf(linspace(0,1,ntrials),[0.05 3 0.4])*.1;   
end
%==========================================================================
% ratios on WEX and WIX
function [ratios1,ratios2]=external_connections_ratios_example3(ntrials)
   ratios1=ones(1,ntrials);
   ratios2=gaussmf(linspace(0,1,ntrials),[0.1 0.8])*.45;
end
%==========================================================================        
function [input,times]=gen_input()

   dt              = 1e-3;    % 1ms
   times           = dt:dt:7; % 7 sec 
   amplitude       = 1.5;     % 1.5 pulses per second
   onofftime       = [3000 3500];

   code.type=4; % piecewise
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

function plot_FRF(allmu,times,ratios1,ratios2,clim)
    % envelope
    allmu_env=allmu;
    for i=1:length(ratios1)
      [pks,locs] = findpeaks(allmu(i,:));
      allmu_env(i,:) = interp1(locs,pks,times*1000);
    end

    figure;
    imagesc(times-3,1:length(ratios1),allmu_env-mean(allmu_env(1,1000:2000)));
    set(gca,'ydir','reverse');
    hold on;
    plot([1;1]*[3 3.5]-3,ylim,'r')
    ylim([0.5 length(ratios1)+0.5])
    xlim([2.5 4]-3)
    set(gca,'xtick',[0 0.5 1],'xticklabel',[0 0.5 1])
    xlabel('Time(s)')
    colorbar
    colormap('jet');caxis(clim);
    title('Simulated onset- and offset-FRFs')

    % plot ratios
    pos=get(gca,'position');
    set(gca,'ytick',[],'position',pos+[0.02 0 0 0 ])
    h=axes('position',[pos(1)-0.08 pos(2) 0.09 pos(4)]);
    hold on
    bar(ratios1,1,'facecolor',[255 102 0]/255,'EdgeColor',[255 102 0]/255);
    bar(ratios2,1,'facecolor','g','EdgeColor','g');
    xlim([0.5 length(ratios1)+0.5])
    view([-90 -90])
    box on; set(gca,'layer','top')
    set(gca,'xtick',[],'ytick',[0 1],'yticklabel',{'0' '100%'})
    ylabel('Input ratio')
end

