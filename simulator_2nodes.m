function [mu,mv,u,v]=simulator_2nodes(input,g0)

N=2; % number of nodes
g=initialize_g(N);


if isfield(g0,'ratio1') && isfield(g0,'ratio2') 
    g.ratio1=g0.ratio1; 
    g.ratio2=g0.ratio2; 
end
if isfield(g0,'condition') 
    g.condition=g0.condition; 
end
if isfield(g0,'c') 
    g.c=g0.c; 
end
p=initialize_p(g);

[WEE,WIE,WEI,WII,WEX,WIX,a]=gen_network(N,g);

    %-------------------start simulation----------------------------- 
    %  u: PSP of E population
    %  v: PSP of I population
    %  mu: firing rate of E population
    %  mv: firing rate of I population
    %  g.state_ue1: PSP of E population contributed by E population (E->E)
    %  g.state_ui1: PSP of E population contributed by I population (I->E)
    %  g.state_ve1: PSP of I population contributed by E population (E->I)
    %  g.state_vi1: PSP of I population contributed by E population (I->I)
    %  W: [N(to) x N(from)]
    %-------------------data initialization-------------------------- 
    u=zeros(N,length(input)); % membrame potential of E populations
    v=zeros(N,length(input)); % membrame potential of I populations
    mu=zeros(N,length(input)); % membrame potential of E populations
    mv=zeros(N,length(input)); % membrame potential of I populations
    
    for t=1:length(input)-1
        if mod(t,1000)==0
            disp([num2str(t) '/' num2str(length(input))])
        end
        u(:,t+1)=g.state_ue1-g.state_ui1;
        v(:,t+1)=g.state_ve1-g.state_vi1;
        mu(:,t+1)=sigm(u(:,t+1),g.v0,g.e0,g.r);% firing rate of E population
        mv(:,t+1)=sigm(v(:,t+1),g.v0,g.e0,g.r);% firing rate of I population

        % E -> E
        g.state_due1 = g.state_ue2;
        g.state_due2 = (g.He./g.taue).*(g.cc*(a.*WEE)*mu(:,t+1)+ WEX*input(:,t) + g.B) - (2./g.taue).*g.state_ue2 -(1./g.taue.^2).*g.state_ue1;  
        % I -> E
        g.state_dui1 = g.state_ui2;
        g.state_dui2 = (g.Hi./g.taui).*(g.cc*WEI*mv(:,t+1)) - (2./g.taui).*g.state_ui2 - (1./g.taui.^2).*g.state_ui1;
        % E -> I
        g.state_dve1 = g.state_ve2;
        g.state_dve2 = (g.He./g.taue).*(g.cc*WIE*mu(:,t+1) + WIX*input(:,t)) - (2./g.taue).*g.state_ve2 - (1./g.taue.^2).*g.state_ve1;
        % I -> I
        g.state_dvi1 = g.state_vi2;
        g.state_dvi2 = (g.Hi./g.taui).*(g.cc*WII*mv(:,t+1)) - (2./g.taui).*g.state_vi2 - (1./g.taui.^2).*g.state_vi1;


        g.state_ue1=g.state_ue1+g.state_due1*g.dt;
        g.state_ue2=g.state_ue2+g.state_due2*g.dt;
        g.state_ui1=g.state_ui1+g.state_dui1*g.dt;
        g.state_ui2=g.state_ui2+g.state_dui2*g.dt;
        g.state_ve1=g.state_ve1+g.state_dve1*g.dt;
        g.state_ve2=g.state_ve2+g.state_dve2*g.dt;
        g.state_vi1=g.state_vi1+g.state_dvi1*g.dt;
        g.state_vi2=g.state_vi2+g.state_dvi2*g.dt;    
        %% plasticity (WEE adaptation) >>>  
        if p.adaptation.on
          da=(1-a)/p.adaptation.taua - p.adaptation.k*a.*repmat(mu(:,t+1)',N,1);
          a=a+da; 
        end 

    end 

end


%==========================================================================  
function g=initialize_g(N,g)

    if nargin<2 % if no g
        % Sigmoid function
        g.v0    = 6e-3; % 6mV
        g.e0    = 2.5;  % 5 s^(-1)
        g.r     = 560;  % 0.56 mV^(-1)

        % Membrane time constant 
        g.taue  = 10*ones(N,1)*1e-3;   % 10 ms   
        g.taui  = g.taue*2;            % 20 ms

        % Averaging synaptic gain
        g.He = 3.25*10*1e-6./g.taue; %to keep He*taue = constant 3.25*10 
        g.Hi = 22*20*1e-6./g.taui;   %to keep Hi*taui = constant 22*20;

        % Averaging number of synaptic constants 
        g.cc= 135;
        g.nNodes=N;

        % Background input (pps)
        g.B=(220/5)*2.5; 

        % Time resolution
        g.dt     = 1e-3;
        
        % Default external connections
        g.ratio1=1; 
        g.ratio2=0; 
        
        % Default condition
        g.condition=1;
        
        % Default inter-column connections
        g.c=[0 0 0 0 0 0 0 0];
    end

    g.state_ue1 =zeros(g.nNodes,1);  % [Nx1] 
    g.state_due1 =zeros(g.nNodes,1); 
    g.state_ue2 =zeros(g.nNodes,1);  
    g.state_due2 =zeros(g.nNodes,1); 

    g.state_ui1 =zeros(g.nNodes,1);  
    g.state_dui1 =zeros(g.nNodes,1); 
    g.state_ui2 =zeros(g.nNodes,1);  
    g.state_dui2 =zeros(g.nNodes,1); 

    g.state_ve1 =zeros(g.nNodes,1);  
    g.state_dve1 =zeros(g.nNodes,1); 
    g.state_ve2 =zeros(g.nNodes,1);  
    g.state_dve2 =zeros(g.nNodes,1); 

    g.state_vi1 =zeros(g.nNodes,1);  
    g.state_dvi1 =zeros(g.nNodes,1); 
    g.state_vi2 =zeros(g.nNodes,1);  
    g.state_dvi2 =zeros(g.nNodes,1); 
end
%==========================================================================  
function p=initialize_p(g)
    % WEE adaptation
    if g.condition==4
      p.adaptation.on = 1;  % 1:ON, 0:OFF 
    else
      p.adaptation.on = 0;     
    end
    p.adaptation.taua=200;  %1.6 sec
    p.adaptation.k=2*1e-3; % constant for adaptation
                            % controlling the falling speed
                            % controlling the min value of a(t) during adaptation

    if p.adaptation.on, disp('WEE adaptation ON'); else, disp('WEE adaptation OFF'); end
end
%==========================================================================  
function [WEE,WIE,WEI,WII,WEX,WIX,a]=gen_network(N,g)
% input: 
%   N: number of nodes
%  
% output:
%   WEE [NxN]                  ****no negative values***
%   WIE [NxN] to I from E      ****no negative values***
%   WEI [NxN] to E from I      ****no negative values***
%   WII [NxN]                  ****no negative values***
%

ifirstLayer=1; % index of first-layer nodes
isecondLayer=2; % index of second-layer nodes
%----connection strength (node)----
%---- good ----
cWEE=0.8;
cWIE=0.6;
cWEI=0.2;           
cWII=0.05;    


%----connection strength (layer 1 to layer 2)----
cWE2E1=g.c(1)/length(ifirstLayer);
cWI2E1=g.c(2)/length(ifirstLayer);
cWE2I1=g.c(3)/length(ifirstLayer);
cWI2I1=g.c(4)/length(ifirstLayer);

%----connection strength (layer 2 to layer 1)----
cWE1E2=g.c(5)/length(ifirstLayer);
cWI1E2=g.c(6)/length(ifirstLayer);
cWE1I2=g.c(7)/length(ifirstLayer);
cWI1I2=g.c(8)/length(ifirstLayer);


%------------WEE------------[NxN] to E from E 
WEE=eye(N)*cWEE;
WEE(isecondLayer,ifirstLayer)=cWE2E1; 
WEE(ifirstLayer,isecondLayer)=cWE1E2; 
a=ones(N); % adaptation term
%------------WEI------------[NxN] to I from E 
WIE=eye(N)*cWIE;
WIE(isecondLayer,ifirstLayer)=cWI2E1; 
WIE(ifirstLayer,isecondLayer)=cWI1E2; 
%------------WIE------------[NxN] to E from I 
WEI=eye(N)*cWEI;
WEI(isecondLayer,ifirstLayer)=cWE2I1; 
WEI(ifirstLayer,isecondLayer)=cWE1I2; 
%------------WII------------[NxN] to I from I 
WII=eye(N)*cWII;
WII(isecondLayer,ifirstLayer)=cWI2I1; 
WII(ifirstLayer,isecondLayer)=cWI1I2; 
%------------WEX------------[Nx1] to E from input 
WEX=(220/5)*[g.ratio1;g.ratio2];      
%------------WIX------------[Nx1] to I from input 
WIX=WEX*0.5;  

end

%==========================================================================  
% v: potential
% m: firing rate
% v0    = 6e-3;
% e0    = 2.5;
% r     = 560;
function m=sigm(v,v0,e0,r)
  m = (2*e0) ./ (1+exp(r*(v0-v)));
end