% To reproduce MMN in roving paradigm
% see /media/vchien/V2/_3.meeting/20180314_onoff_3nodes_rovingParadigm/20180314_onoff_3nodes.odp
%
% 2_1 network
%   display 6x6 array of node3 responses (input1 strength 0:0.2:1  x input2 strength 0:0.2:1)
% layer1:
%    input 1----> node1
%    input 2----> node2
% layer2:
%    (no input)-> node3
%
%    input1----> node1 <---->        (type2: sustained on)
%                             node3
%    input2----> node2 <---->        (type6: inhibited off)
% 
% see 6 types of on-off responses at node3
function Fig5()
    % Input: prolonged stimulus
    g.onofftime       = [2 4 6.5 9 11.5 13.5]*1000; % to match paper '(2016). Brain responses in humans reveal ideal observer-like sensitivity...'
    [input,times]=gen_input(g);
    %---------------Fig. 6a-------------------
    % Inter-column connections
    % 
    g.c31=[1 4 1 2]*0.1; %(node 1->3)      sustained-On&Off
    g.c13=[1 1 2 1]*0.1; %(node 3->1)          
    g.c32=[1 4 1 2]*0.1; %(node 2->3)      sustained-On&Off
    g.c23=[1 1 2 1]*0.1; %(node 3->2)            
    g.c21=[4 4 2   2  ]*0.1; %(node 1->2)      inter-feature
    g.c12=[0 0 2.5 2.5]*0.1; %(node 2->1)  
    
    g.condition=1; 
    % Run simulation
    [mu,mv,u,v,R]=simulator_3nodes(input,g);
    plot_result(R,times,g);
    plot_result2(times,3,input,mu,mv,u,v,g,g.onofftime,R)
end

%========================================================================== 
% generating roving paradigm
function [input,times]=gen_input(g)

   dt              = 1e-3;     % 1ms
   times           = dt:dt:16; % 16 sec 
   amplitude       = 1.5;      % 1.5 pulses per second
   input=zeros(2,length(times)); 
   
   
   code.type=4; % piecewise
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
%==========================================================================
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
    %xlabel('Time (ms)');
    ylabel('Simulated MEG');
    legend('Envelope','Simulated MEG')
    subplot(212);
    onofftime=g.onofftime-9000;
    plot(times*1000-9000,R_env,'color',[255 145 0]/255,'linewidth',3); hold on
    plot(times*1000-9000,R,'k'); 
    plot([1;1]*onofftime(4),ylim,'k:'); xlim([-500 5000]);
    plot([1;1]*onofftime(5),ylim,'k:'); xlim([-500 5000]);
    plot([1;1]*onofftime(6),ylim,'k:'); xlim([-500 5000]);
    title('REG-RAND'); 
    xlabel('Time (ms)');ylabel('Simulated MEG');
    %legend('Envelope','Simulated MEG')
    set(gcf,'position',[0 0 600 600])
end
%==========================================================================
function plot_result2(times,N,input,mu,mv,u,v,g,onofftime,R)

       figure;
%         subplot(4,1,1);hold on 
%         for i=1:N
%             plot(times,inputE(i,:)+5*(i-1),'k','LineWidth',2)
%             plot(times,inputI(i,:)+5*(i-1),'r')
%         end
%         plot([1;1]*onofftime*0.001,ylim,'b'); % blue lines
%         plot(xlim,[1;1]*[0:N-1]*5,'k'); % hori lines 
% 
%         ylabel('firing rate')
%         set(gca,'ytick',[0 4 5 9 10 14],'yticklabel',[0 4 0 4 0 4])
%         
%         %title({sprintf('node1:[%g,%g,%g,%g]',g.c(1:4)),sprintf('node2:[%g,%g,%g,%g]',g.c(5:8)), sprintf('tau_e_,_i=[%g,%g]ms', g.taue(1)*1000,g.taui(1)*1000),'Inputs (to node 1,2)'})
%         ylim([0 15])
        subplot(4,1,1); %plot_network(g)

        subplot(4,1,2);hold on
        for i=1:N
            plot(times,mu(i,:)'+5*(i-1),'k','LineWidth',2);
            plot(times,mv(i,:)'+5*(i-1),'r','LineWidth',1);    
        end
        plot([1;1]*onofftime*0.001,ylim,'b'); % blue lines
        plot(xlim,[1;1]*[0:N-1]*5,'k'); % hori lines 
        set(gca,'ytick',[0 4 5 9 10 14],'yticklabel',[0 4 0 4 0 4])
        ylabel('firing rate')
        legend('E','I')
        title('Responses (node 1,2,3)')
        ylim([0 15])

        subplot(4,1,3);hold on
        scale=max([u(:) ;v(:)])-min([u(:) ;v(:)]);
        for i=1:N
            plot(times,u(i,:)'+scale*(i-1),'k','LineWidth',2);
            plot(times,v(i,:)'+scale*(i-1),'r','LineWidth',1);        
        end
        plot([1;1]*onofftime*0.001,ylim,'b'); % blue lines
        plot(xlim,[1;1]*[0:N-1]*scale,'k'); % hori lines
        
        set(gca,'ytick',[0 0.01 scale scale*2],'yticklabel',[0 0.01 0 0])
        ylabel('potential')
        % legend('I','E')
        title('Responses (node 1,2,3)')
        ylim([min([u(:) ;v(:)]), scale*2+max([u(:) ;v(:)])])
        xlabel('time (sec)')
        
        subplot(4,1,4);hold on;
        ylimit=ceil(max(R)*1.1);
        area(times,(input(1,:)~=0)*ylimit,'edgecolor',[255 102 0]/255,'facecolor',[255 102 0]/255)
        area(times,(input(2,:)~=0)*ylimit,'edgecolor',[0 153 0]/255,'facecolor',[0 153 0]/255)
        text(2,ylimit,'tone 1','HorizontalAlignment','center','VerticalAlignment','top'); text(9,ylimit,'tone 1','HorizontalAlignment','center','VerticalAlignment','top');
        text(4,ylimit,'tone 2','HorizontalAlignment','center','VerticalAlignment','top'); text(7,ylimit,'tone 2','HorizontalAlignment','center','VerticalAlignment','top');
        plot(times,R,'k','LineWidth',2); 
        plot([1;1]*onofftime*0.001,ylim,'b'); % blue lines
        ylim([0 ylimit])
        if g.soa==0
          title('simulated MEG')
        else
          title(sprintf('simulated MEG (SOA:%g,duration:%g ms)',g.soa,g.duration));
        end
        xlabel('time (sec)')
        set(gcf,'position',[0 0 900 960]);


end
