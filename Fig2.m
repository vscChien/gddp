% Shih-Cheng Chien, Burkhard Maess, and Thomas R. Knösche. "A generic deviance detection principle for cortical On/Off responses, omission response, and mismatch negativity." BioRxiv (2019): 582437.
function Fig2()
  %-------------load data--------------
  disp('Load W solutions (condition 1)')
  load('solutions_condition1(downsampled).mat');
  fprintf('Inc-None   : %d solutions\n',sum(groups1==1 & groups2==0));
  fprintf('Inc-On     : %d solutions\n',sum(groups1==1 & groups2==1));
  fprintf('Inc-Off    : %d solutions\n',sum(groups1==1 & groups2==2));
  fprintf('Inc-OnOff  : %d solutions\n',sum(groups1==1 & groups2==3));
  fprintf('Dec-None   : %d solutions\n',sum(groups1==2 & groups2==0));
  fprintf('Dec-On     : %d solutions\n',sum(groups1==2 & groups2==1));
  fprintf('Dec-Off    : %d solutions\n',sum(groups1==2 & groups2==2));
  fprintf('Dec-OnOff  : %d solutions\n',sum(groups1==2 & groups2==3));
  
  disp('Run t-SNE...(around 1 minute)')
  rng default 
  Y = tsne(allsolutions,'Algorithm','barneshut','Distance','euclidean'); 


  % eight exemplary types in Fig.2cd
  iloc=[896;  % Inc-None
        1220; % Inc-On
        3013; % Inc-Off
        3219; % Inc-OnOff
        3577; % Dec-None
        4280; % Dec-On
        5819; % Dec-Off
        5941];% Dec-OnOff
   
  plot_2D(Y,groups1,groups2,iloc);
  plot_responses(allsolutions,iloc);
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
function plot_2D(Y,groups1,groups2,iloc)

 figure; 
 gscatter(Y(:,1),Y(:,2),groups1,[0.9 0.5 0.5;...        % Inc
                                 0.5 0.5 0.9],'.',20);  % Dec
 hold on;                            
 scatter(Y(groups2==1,1),Y(groups2==1,2),20,'r.');      % On
 scatter(Y(groups2==2,1),Y(groups2==2,2),20,'g.');      % Off
 scatter(Y(groups2==3,1),Y(groups2==3,2),20,'b.');      % OnOff
 rectangle('Position',[-40 17 15 20]);
 title('Distribution of On/Off types in W space'); 
 legend('Inc','Dec','On','Off','OnOff'); 
 xlim([-1 1]*70);ylim([-1 1]*70);
 set(gca,'ytick',[],'xtick',[])
 box on
 set(gcf,'position',[0 0 600 600])
 
 figure; 
 gscatter(Y(:,1),Y(:,2),groups1,[0.9 0.5 0.5;...        % Inc
                                 0.5 0.5 0.9],'.',20);  % Dec   
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
 set(gca,'ytick',[],'xtick',[]);
 box on
 set(gcf,'position',[0 0 200 200]);

end
%--------------------------------------------------------------------------
function plot_responses(allsolutions,iloc)
 T={'Inc-None','Inc-On','Inc-Off','Inc-OnOff','Inc-None','Inc-On','Inc-Off','Inc-OnOff'};
 signals=zeros(length(iloc),7000);
 signals_env=zeros(length(iloc),7000);
 
 [input,times]=gen_input(); % prolonged stimulus (2 sec)
 for i=1:length(iloc)
    g.c=allsolutions(iloc(i),:);
    [mu,~,~,~]=simulator_2nodes(input,g);
    signals(i,:)=mu(2,:);
    [pks,locs] = findpeaks(signals(i,:));
    signals_env(i,:) = interp1(locs,pks,times*1000); % envelope
 end
 figure;
 for i=1:length(iloc)
     subplot(8,1,i);hold on
     area([0 3 3 5 5 6]',[0 0 5 5 0 0]','facecolor',[1 1 1]*0.5,'edgecolor',[1 1 1]*0.5)
     plot(times,signals_env(i,:),'k','linewidth',2);
     plot(times,signals(i,:),'k');
     xlim([2 6]); ylim([0 5]);
     title(['(' num2str(i) ') ' T{i}]);
     axis off
 end
 
 set(gcf,'position',[0 0 200 900]);

end
