function Fig7()

    % Input: prolonged stimulus
    [input,times]=gen_input();
    %---------------Fig. 7a-------------------
    % Inter-column connections
    g.c=[0 1 1 2 1 1 2 1]*0.1;   % inhibited-Off
    g.condition=1; 
    % Run simulation
    [mu,mv,u,v]=simulator_2nodes(input,g);
    plot_result(mu,mv,u,v,times);
    %---------------Fig. 7b-------------------
    % Inter-column connections
    g.c=[1 1 1 2 2 1 2 2]*0.1;   % sustained-Off
    g.condition=1; 
    % Run simulation
    [mu,mv,u,v]=simulator_2nodes(input,g);
    plot_result(mu,mv,u,v,times);

end


%========================================================================== 
function plot_result(mu,mv,u,v,times)

figure;
subplot(2,1,1);hold on
area([3 3 5 5],[0 10 10 0],'facecolor',[1 1 1]*0.7,'linestyle','none');%green color stimulus
h1=plot(times,mu(1,:)','k','LineWidth',2);
h2=plot(times,mv(1,:)','r','LineWidth',1);

plot(times,mu(2,:)'+5,'k','LineWidth',2);
plot(times,mv(2,:)'+5,'r','LineWidth',1);

xlim([1 7]);ylim([0 10]);
plot(xlim,[1;1]*[0,5],'k'); % black lines 
set(gca,'ytick',[0 1 5 6],'yticklabel',[0 1 0 1])
%set(gca,'ytick',[],'yticklabel',[])
%ylabel('Col.1  Col.2')

legend([h1,h2],{'E pop.','I pop.'})
text(1.1,4,'column 1')
text(1.1,9,'column 2')
title('Column responses')
%xlabel('Time (sec)')
ylabel('Firing rate (pps)')
set(gca,'layer','top')
box on


subplot(2,1,2);hold on
area([3 3 5 5],[0 5 5 0],'facecolor',[1 1 1]*0.7,'linestyle','none');%green color stimulus
h1=plot(times,u(1,:)'+0.02,'k','LineWidth',2);
h2=plot(times,v(1,:)'+0.02,'r','LineWidth',1);
plot(times,u(2,:)'+0.05,'k','LineWidth',2);
plot(times,v(2,:)'+0.05,'r','LineWidth',1);
plot(xlim,[1;1]*[0.02,0.05],'k'); % blue lines
xlim([1 7]);ylim([0 0.07]);

set(gca,'ytick',[0 0.01 0.03 0.04]+0.02,'yticklabel',[0 0.01 0 0.01])
%legend([h1,h2],{'E pop.','I pop.'})
text(1.1,0.035,'column 1')
text(1.1,0.065,'column 2')
%title('Column responses')
xlabel('Time (sec)')
ylabel('PSP (mV)')
set(gca,'layer','top')
box on

set(gcf,'position',[0 0 460 520])

end

%==========================================================================        
function [input,times]=gen_input()

   dt              = 1e-3;    % 1ms
   times           = dt:dt:7; % 7 sec 
   amplitude       = 1.5;     % 1.5 pulses per second
   onofftime       = [3000 5000];

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