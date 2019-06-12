% Shih-Cheng Chien, Burkhard Maess, and Thomas R. Knösche. "A generic deviance detection principle for cortical On/Off responses, omission response, and mismatch negativity." BioRxiv (2019): 582437.
function Fig5()
    % input: prolonged stimulus
    g.onofftime       = [2 4 6.5 9 11.5 13.5]*1000; 
    [input,times]=gen_input(g);
    %---------------Fig. 6a-------------------
    % inter-column connections
    g.c31=[1 4 1 2]*0.1; %(node 1->3)      Inc-OnOff
    g.c13=[1 1 2 1]*0.1; %(node 3->1)          
    g.c32=[1 4 1 2]*0.1; %(node 2->3)      Inc-OnOff
    g.c23=[1 1 2 1]*0.1; %(node 3->2)            
    g.c21=[4 4 2   2  ]*0.1; %(node 1->2)      
    g.c12=[0 0 2.5 2.5]*0.1; %(node 2->1)  
    
    g.condition=1; 
    % run simulation
    [~,~,~,~,R]=simulator_3nodes(input,g);
    plot_result(R,times,g);
end
%--------------------------------------------------------------------------  
function [input,times]=gen_input(g)

   dt              = 1e-3;     % 1 msec
   times           = dt:dt:16; % 16 sec 
   amplitude       = 1.5;      % 1.5 pps
   input=zeros(2,length(times)); 

   code.n=1;
   code.tlength=length(times);
   code.dt=dt; 
   code.peak=amplitude; 
   code.soa=20000;
   code.duration=diff(g.onofftime);
   code.risetime=10;
   code.figureon=0;
   
   %RAND
   code.offset=g.onofftime(1);                             
   code.duration=g.onofftime(2)-g.onofftime(1)+code.risetime;
   input1_1=gen_input_core(code)'; 
   %REG
   code.offset=g.onofftime(2);                             
   code.duration=g.onofftime(3)-g.onofftime(2)+code.risetime;
   input2_1=gen_input_core(code)'; 
   %REG
   code.offset=g.onofftime(4);                             
   code.duration=g.onofftime(5)-g.onofftime(4)+code.risetime;
   input2_2=gen_input_core(code)'; 
   %RAND
   code.offset=g.onofftime(5);                            
   code.duration=g.onofftime(6)-g.onofftime(5)+code.risetime;
   input1_2=gen_input_core(code)'; 
   
   
   input(1,:)=input1_1+input1_2;%---<RAND>-----------<RAND>---
   input(2,:)=input2_1+input2_2;%--------<REG>---<REG>--------

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
function plot_result(R,times,g)
    
    [pks,locs] = findpeaks(R);
    R_env = interp1(locs,pks,times*1000);
    figure;
    subplot(211);
    onofftime=g.onofftime-2000;
    plot(times*1000-2000,R_env,'color',[255 145 0]/255,'linewidth',3); hold on
    plot(times*1000-2000,R,'k'); 
    plot([1;1]*onofftime(1),ylim,'k:'); xlim([-500 5000]);
    plot([1;1]*onofftime(2),ylim,'k:'); xlim([-500 5000]);
    plot([1;1]*onofftime(3),ylim,'k:'); xlim([-500 5000]);
    title('RAND-REG'); 
    ylabel('Simulated MEG');
    legend('Envelope','Simulated MEG');
    subplot(212);
    onofftime=g.onofftime-9000;
    plot(times*1000-9000,R_env,'color',[255 145 0]/255,'linewidth',3); hold on
    plot(times*1000-9000,R,'k'); 
    plot([1;1]*onofftime(4),ylim,'k:'); xlim([-500 5000]);
    plot([1;1]*onofftime(5),ylim,'k:'); xlim([-500 5000]);
    plot([1;1]*onofftime(6),ylim,'k:'); xlim([-500 5000]);
    title('REG-RAND'); 
    xlabel('Time (ms)');ylabel('Simulated MEG');
    set(gcf,'position',[0 0 600 600]);
end

