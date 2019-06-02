function Fig2()

  %------load data---------------------
  disp('Load W solutions (condition 1)')
  load('solutions_condition1(downsampled).mat');
  fprintf('Sustained-Pure   : %d solutions\n',sum(groups1==1 & groups2==0));
  fprintf('Sustained-On     : %d solutions\n',sum(groups1==1 & groups2==1));
  fprintf('Sustained-Off    : %d solutions\n',sum(groups1==1 & groups2==2));
  fprintf('Sustained-On&Off : %d solutions\n',sum(groups1==1 & groups2==3));
  fprintf('Inhibited-Pure   : %d solutions\n',sum(groups1==2 & groups2==0));
  fprintf('Inhibited-On     : %d solutions\n',sum(groups1==2 & groups2==1));
  fprintf('Inhibited-Off    : %d solutions\n',sum(groups1==2 & groups2==2));
  fprintf('Inhibited-On&Off : %d solutions\n',sum(groups1==2 & groups2==3));
  
  disp('Run t-SNE...(around 1 minute)')
  rng default % for reproducibility 
  Y = tsne(allsolutions,'Algorithm','barneshut','Distance','euclidean'); 


  iloc=[896;  % sustained pure
        1220; % sustained On
        3013; % sustained Off
        3219; % sustained On&Off
        3577; % inhibited pure
        4280; % inhibited On
        5819; % inhibited Off
        5941];% inhibited On&Off
   
  plot_2D(Y,groups1,groups2,iloc);
  plot_responses(allsolutions,iloc);
end
 
 

%--------------------------------------------------------------------------
function plot_2D(Y,groups1,groups2,iloc)

 figure; 
 gscatter(Y(:,1),Y(:,2),groups1,[0.9 0.5 0.5;...        % sustained
                                 0.5 0.5 0.9],'.',20);  % inhibited
 hold on;                            
 scatter(Y(groups2==1,1),Y(groups2==1,2),20,'r.');      % On
 scatter(Y(groups2==2,1),Y(groups2==2,2),20,'g.');      % Off
 scatter(Y(groups2==3,1),Y(groups2==3,2),20,'b.');      % On&Off
 rectangle('Position',[-40 17 15 20]);
 title('Distribution of On/Off types in W space'); 
 legend('Sustained','Inhibited','On','Off','On&Off'); 
 xlim([-1 1]*70);ylim([-1 1]*70);
 set(gca,'ytick',[],'xtick',[])
 box on
 set(gcf,'position',[0 0 600 600])
 
 figure; 
 gscatter(Y(:,1),Y(:,2),groups1,[0.9 0.5 0.5;...        % sustained
                                 0.5 0.5 0.9],'.',20);  % inhibited   
 hold on;                            
 scatter(Y(groups2==1,1),Y(groups2==1,2),20,'r.');      % On
 scatter(Y(groups2==2,1),Y(groups2==2,2),20,'g.');      % Off
 scatter(Y(groups2==3,1),Y(groups2==3,2),20,'b.');      % OnOff
 for i=1:length(iloc)
    scatter(Y(iloc(i),1),Y(iloc(i),2),'ok');
    text(Y(iloc(i),1),Y(iloc(i),2),num2str(i),'VerticalAlignment','bottom','HorizontalAlignment','center');   
 end
 title('zoom in'); 
 legend off
 xlim([-40 -25]);ylim([17 37]);
 set(gca,'ytick',[],'xtick',[])
 box on
 set(gcf,'position',[0 0 200 200])

end

%--------------------------------------------------------------------------
function plot_responses(allsolutions,iloc)
 T={'Sustained Pure','Sustained On','Sustained Off','Sustained On&Off','Inhibited Pure','Inhibited On','Inhibited Off','Inhibited On&Off'};
 signals=zeros(length(iloc),7000);
 signals_env=zeros(length(iloc),7000);
 
 % Input: prolonged stimulus
 [input,times]=gen_input(); 
 for i=1:length(iloc)
    g.c=allsolutions(iloc(i),:);
    [mu,~,~,~]=simulator_2nodes(input,g);
    signals(i,:)=mu(2,:);
    [pks,locs] = findpeaks(signals(i,:));
    signals_env(i,:) = interp1(locs,pks,times*1000);
 end
 figure;
 for i=1:length(iloc)
     subplot(8,1,i);hold on
     area([0 3 3 5 5 6]',[0 0 5 5 0 0]','facecolor',[1 1 1]*0.5,'edgecolor',[1 1 1]*0.5)
     plot(times,signals_env(i,:),'k','linewidth',2);
     plot(times,signals(i,:),'k');
     xlim([2 6])
     ylim([0 5])
     title(['(' num2str(i) ') ' T{i}])
     axis off
 end
 
 set(gcf,'position',[0 0 200 900])

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
 