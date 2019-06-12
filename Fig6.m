% Shih-Cheng Chien, Burkhard Maess, and Thomas R. Knösche. "A generic deviance detection principle for cortical On/Off responses, omission response, and mismatch negativity." BioRxiv (2019): 582437.
function Fig6()
    % input: prolonged stimulus
    [input,times]=gen_input();
    %---------------Fig. 6a-------------------
    % inter-column connections
    g.c=[0 2 1 2 1 2 2 1]*0.1;   % Dec-OnOff
    g.condition=1; 
    % run simulation
    [mu,mv,u,v]=simulator_2nodes(input,g);
    plot_result(mu,mv,u,v,times);
    %---------------Fig. 6b-------------------
    % inter-column connections
    g.c=[1 4 1 2 1 1 2 1]*0.1;   % Inc-OnOff
    g.condition=1; 
    % run simulation
    [mu,mv,u,v]=simulator_2nodes(input,g);
    plot_result(mu,mv,u,v,times);

end
%--------------------------------------------------------------------------         
function [input,times]=gen_input()

   dt              = 1e-3;        % 1 msec
   times           = dt:dt:7;     % 7 sec 
   amplitude       = 1.5;         % 1.5 pps
   onofftime       = [3000 5000]; % 2 sec

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
function plot_result(mu,mv,u,v,times)

figure;
subplot(2,1,1);hold on
area([3 3 5 5],[0 10 10 0],'facecolor',[1 1 1]*0.7,'linestyle','none');
h1=plot(times,mu(1,:)','k','LineWidth',2);
h2=plot(times,mv(1,:)','r','LineWidth',1);

plot(times,mu(2,:)'+5,'k','LineWidth',2);
plot(times,mv(2,:)'+5,'r','LineWidth',1);

xlim([1 7]);ylim([0 10]);
plot(xlim,[1;1]*[0,5],'k'); 
set(gca,'ytick',[0 1 5 6],'yticklabel',[0 1 0 1]);

legend([h1,h2],{'E pop.','I pop.'});
text(1.1,4,'node 1');
text(1.1,9,'node 2');
title('Population responses');
ylabel('Firing rate (spikes/s)');
set(gca,'layer','top');
box on


subplot(2,1,2);hold on
area([3 3 5 5],[0 5 5 0],'facecolor',[1 1 1]*0.7,'linestyle','none');
h1=plot(times,u(1,:)'+0.02,'k','LineWidth',2);
h2=plot(times,v(1,:)'+0.02,'r','LineWidth',1);
plot(times,u(2,:)'+0.05,'k','LineWidth',2);
plot(times,v(2,:)'+0.05,'r','LineWidth',1);
plot(xlim,[1;1]*[0.02,0.05],'k'); 
xlim([1 7]);ylim([0 0.07]);

set(gca,'ytick',[0 0.01 0.03 0.04]+0.02,'yticklabel',[0 0.01 0 0.01]);
text(1.1,0.035,'node 1');
text(1.1,0.065,'node 2');
xlabel('Time (sec)');
ylabel('PSP (mV)');
set(gca,'layer','top');
box on
set(gcf,'position',[0 0 460 520]);

end
