function Fig8()
    % 4 conditions: 1: adaptation OFF, Wix=1
    %               2: adaptation OFF, Wix=0
    %               3: NMDAblocked: WEE=WEE*0.75, WEI=WEI*0.5
    %               4: adaptationON, Wix=1
    allsolutions_4conds={}; %4x1 cells
    groupsSI_4conds={};     %4x1 cells
    groupsOnOff_4conds={};  %4x1 cells
    onoff_4conds={}; %4x1 cells
    type_4conds={}; %4x1 cells

    for i=1:4 % the 4 conditions
      fprintf('Load W solutions (condition %d)\n',i);
      load(sprintf('solutions_condition%d.mat',i));

      solutions_sustainedPure=solutions(solutions_type==1,:);
      solutions_sustainedOn=solutions(solutions_type==2,:);
      solutions_sustainedOff=solutions(solutions_type==3,:);
      solutions_sustainedOnOff=solutions(solutions_type==4,:);
      solutions_inhibitedPure=solutions(solutions_type==5,:);
      solutions_inhibitedOn=solutions(solutions_type==6,:);
      solutions_inhibitedOff=solutions(solutions_type==7,:);
      solutions_inhibitedOnOff=solutions(solutions_type==8,:);
      type_4conds{i}=solutions_type;

      onoff_4conds{i}=[solutions_sustainedOn;...
                       solutions_sustainedOff;...
                       solutions_sustainedOnOff;...
                       solutions_inhibitedOn;...
                       solutions_inhibitedOff;...
                       solutions_inhibitedOnOff];  


      allsolutions_4conds{i}=[solutions_sustainedPure;...
                              solutions_sustainedOn;...
                              solutions_sustainedOff;...
                              solutions_sustainedOnOff;...
                              solutions_inhibitedPure;...
                              solutions_inhibitedOn;...
                              solutions_inhibitedOff;...
                              solutions_inhibitedOnOff];  
      groupsSI_4conds{i}=[ones(size(solutions_sustainedPure,1),1)*1;...
                          ones(size(solutions_sustainedOn,1),1)*1;...
                          ones(size(solutions_sustainedOff,1),1)*1;...
                          ones(size(solutions_sustainedOnOff,1),1)*1;...
                          ones(size(solutions_inhibitedPure,1),1)*2;...
                          ones(size(solutions_inhibitedOn,1),1)*2;...
                          ones(size(solutions_inhibitedOff,1),1)*2;...
                          ones(size(solutions_inhibitedOnOff,1),1)*2];                   
      groupsOnOff_4conds{i}=[ones(size(solutions_sustainedPure,1),1)*0;...
                             ones(size(solutions_sustainedOn,1),1)*1;...
                             ones(size(solutions_sustainedOff,1),1)*2;...
                             ones(size(solutions_sustainedOnOff,1),1)*3;...
                             ones(size(solutions_inhibitedPure,1),1)*0;...
                             ones(size(solutions_inhibitedOn,1),1)*1;...
                             ones(size(solutions_inhibitedOff,1),1)*2;...
                             ones(size(solutions_inhibitedOnOff,1),1)*3];  

    end
    allsolutions=unique([onoff_4conds{1};onoff_4conds{2};onoff_4conds{3};onoff_4conds{4}],'rows'); %merge the 4 conditions 
    fprintf('%d solutions in total\n',length(allsolutions));

%     disp('Run t-SNE...(around 5 minutes)')
%     rng default % for reproducibility 
%     Y = tsne(allsolutions,'Algorithm','barneshut','Distance','euclidean'); 
% 
%     plot_2D(Y,allsolutions,allsolutions_4conds,groupsSI_4conds,groupsOnOff_4conds);
    plot_table(type_4conds{1},type_4conds{2},'condition I','condition II');
    plot_table(type_4conds{1},type_4conds{3},'condition I','condition III');
    plot_table(type_4conds{1},type_4conds{4},'condition I','condition IV');
    plot_bar(type_4conds);

end

%--------------------------------------------------------------------------
function plot_2D(Y,allsolutions,allsolutions_4conds,groupsSI_4conds,groupsOnOff_4conds)

    for i=1:4 
        [c,idx] = ismember(allsolutions,allsolutions_4conds{i},'rows');
        Ytmp=Y(c,:);
        idx(idx==0)=[];
        groupsSI=groupsSI_4conds{i}(idx,:);
        groupsOnOff=groupsOnOff_4conds{i}(idx,:);
        
        figure; 
        scatter(Y(:,1),Y(:,2),30,[1 1 1]*0.8,'filled');  hold on;       % All
        gscatter(Ytmp(:,1),Ytmp(:,2),groupsSI,[0.9 0.5 0.5;...          % sustained
                                               0.5 0.5 0.9],'.',20);    % inhibited                               
        scatter(Ytmp(groupsOnOff==1,1),Ytmp(groupsOnOff==1,2),20,'r.'); % On
        scatter(Ytmp(groupsOnOff==2,1),Ytmp(groupsOnOff==2,2),20,'g.'); % Off
        scatter(Ytmp(groupsOnOff==3,1),Ytmp(groupsOnOff==3,2),20,'b.'); % On&Off     
        switch i
            case 1
               title({'Distribution of On/Off types in W space' ,'(condition I: default)'}); 
            case 2
               title({'Distribution of On/Off types in W space' ,'(condition II: W^I^X=0)'}); 
            case 3
               title({'Distribution of On/Off types in W space' ,'(condition III: NMDA-r ant.)'}); 
            case 4
               title({'Distribution of On/Off types in W space' ,'(condition IV: adaptation)'}); 
        end       
        legend('Others','Sustained','Inhibited','On','Off','On&Off');  
        xlim([-1 1]*70);ylim([-1 1]*70);
        set(gca,'ytick',[],'xtick',[])
        box on
        set(gcf,'position',[0 0 600 600])
    end
end

%--------------------------------------------------------------------------
function plot_table(cond1,cond2,text1,text2)

    % re-order
    temp=[cond1,cond2];
    temp2=temp;
    temp2(temp==2)=1; %sustainedOn
    temp2(temp==6)=2; %inhibitedOn
    temp2(temp==3)=3; %sustainedOff
    temp2(temp==7)=4; %inhibitedOff
    temp2(temp==4)=5; %sustainedOnOff
    temp2(temp==8)=6; %inhibitedOnOff
    temp2(temp==1)=7; %sustainedPure
    temp2(temp==5)=8; %inhibitedPure
    temp2(temp==0)=9; %others
    mat=zeros(9,9);
    for i=1:9
        for j=1:9
            mat(i,j)=sum(ismember(temp2,[i j],'rows'));
        end
    end
    % normalize to cond1
%   nOnOff=sum(temp2(:,1)<=6); % number of On/Off/On&Off
    mat2=round(mat*100/length(temp),2); % in percentage

    figure; imagesc(mat2); 
    hold on;
    plot([1 ;1]*[1.5:8.5],ylim,'k')
    plot(xlim,[1 ;1]*[1.5:8.5],'k')
    rectangle('Position',[6.5 0.5 3 6],'edgecolor','c','linewidth',2);
    rectangle('Position',[0.5 6.5 6 3],'edgecolor','m','linewidth',2);
    for i=1:9
        for j=1:9
            if mat(i,j)~=0
              text(j,i,num2str(mat2(i,j)),'HorizontalAlignment','center','color','r')
            end
        end
    end
    set(gca,'ticklength',[0 0],'xaxisLocation','top'...
        ,'xticklabel',{'On_s','On_i','Off_s','Off_i','OnOff_s','OnOff_i','Pure_s','Pure_i','Others'}...
        ,'yticklabel',{'On_s','On_i','Off_s','Off_i','OnOff_s','OnOff_i','Pure_s','Pure_i','Others'})
    colormap('gray')
    caxis([0 1])
    ylabel(text1)
    xlabel(text2)
end

%--------------------------------------------------------------------------
function plot_bar(type_4conds)

    temp=[type_4conds{1}';...
          type_4conds{2}';...
          type_4conds{3}';...
          type_4conds{4}'];

    temp2=zeros(size(temp));
    temp2(temp==2)=1; %sustained On
    temp2(temp==3)=2; %sustained Off
    temp2(temp==4)=3; %sustained OnOff
    temp2(temp==6)=4; %inhibited On
    temp2(temp==7)=5; %inhibited Off
    temp2(temp==8)=6; %inhibited OnOff

    bartable=zeros(4,6); % [4 condition X 6 types]
    for i=1:4
        for j=1:6
            bartable(i,j)=sum(temp2(i,:)==j);
        end
    end

    % normalize to cond1
    %bartable=bartable/sum(bartable(1,:));
    bartable=bartable*100/length(temp2); %percentage
    figure;
    b=bar(bartable,'grouped');
    xlabels = {'Condition I\newline(Default)','Condition II\newline(W^I^X=0)','Condition III\newline(NMDA-r ant.)','Condition IV\newline(Adaptation)'};
    set(gca,'xticklabel',xlabels,'ytick',[0 1 2 3 4],'YGrid','on');
    ylabel('Proportion (%)')

    b(1).FaceColor=[1 0 0];
    b(1).EdgeColor=[0.9 0.5 0.5];
    b(1).LineWidth=1.5;
    b(2).FaceColor=[0 1 0];
    b(2).EdgeColor=[0.9 0.5 0.5];
    b(2).LineWidth=1.5;
    b(3).FaceColor=[0 0 1];
    b(3).EdgeColor=[0.9 0.5 0.5];
    b(3).LineWidth=1.5;
    b(4).FaceColor=[1 0 0];
    b(4).EdgeColor=[0.5 0.5 0.9];
    b(4).LineWidth=1.5;
    b(5).FaceColor=[0 1 0];
    b(5).EdgeColor=[0.5 0.5 0.9];
    b(5).LineWidth=1.5;
    b(6).FaceColor=[0 0 1];
    b(6).EdgeColor=[0.5 0.5 0.9];
    b(6).LineWidth=1.5;

    legend({'On_s','Off_s','OnOff_s','On_i','Off_i','OnOff_i'},'location','northwest')
    title('On/Off types')
    ylim([0 4])
end

