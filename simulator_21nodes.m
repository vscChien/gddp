function [mu,mv,u,v,R]=simulator_21nodes(input,g0)

N=21; % number of nodes
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
p=initialize_p;

[WEE,WIE,WEI,WII,WEX,WIX,gwE,gwI]=gen_network(N,g);

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
u=zeros(N,length(input));  
v=zeros(N,length(input));  
mu=zeros(N,length(input)); 
mv=zeros(N,length(input)); 

for t=1:length(input)-1
    if mod(t,1000)==0
        disp([num2str(t) '/' num2str(length(input))])
    end
    u(:,t+1)=g.state_ue1-g.state_ui1;
    v(:,t+1)=g.state_ve1-g.state_vi1;
    mu(:,t+1)=sigm(u(:,t+1),g.v0,g.e0,g.r);
    mv(:,t+1)=sigm(v(:,t+1),g.v0,g.e0,g.r);

    % E -> E
    g.state_due1 = g.state_ue2;
    g.state_due2 = (g.He./g.taue).*(g.cc*WEE*mu(:,t+1)+ WEX*input(:,t) + g.B) - (2./g.taue).*g.state_ue2 -(1./g.taue.^2).*g.state_ue1;  
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

    % ----------plasticity ----------
    if p.resonance.on && t>p.resonance.memoryLength 
         %-------fading window-------
         fadingWindow=linspace(0,1,p.resonance.memoryLength);  
         %-------get covariance-------
         tmpu=mu(:,t-p.resonance.memoryLength+1:t);
         tmpu=tmpu-repmat(mean(tmpu,2),[1,size(tmpu,2)]); % remove mean 
         tmpu=tmpu.*repmat(fadingWindow,N,1);
         covMatrixu=cov(tmpu');
         covMatrixu=covMatrixu./max(abs(covMatrixu(:)));
         %-------update WEE & WEI-------
         dWEE=(-WEE + p.resonance.WEE_rate*abs(covMatrixu.*(covMatrixu> 0)).*gwE)/1000; 
         dWEE(eye(N)==1)=0; dWEE(:,N)=0; dWEE(N,:)=0; 
         WEE=max(WEE+dWEE,0);
         dWEI=(-WEI + p.resonance.WEI_rate*abs(covMatrixu.*(covMatrixu< 0)).*gwI)/1000; 
         dWEI(eye(N)==1)=0;dWEI(:,N)=0; dWEI(N,:)=0; 
         WEI=max(WEI+dWEI,0);   
    end
end 
%----------simulated MEG-----------
MEGweight=([ones(1,N-1), (N-1)*3]);
MEGweight=MEGweight/sum(MEGweight);
R=MEGweight*(WEE*mu + WEI*mv);

end


%--------------------------------------------------------------------------  
function g=initialize_g(N,g)

    if nargin<2 % if no g
        % Sigmoid function
        g.v0    = 6e-3; % 6mV
        g.e0    = 2.5;  % 5 s^(-1)
        g.r     = 560;  % 0.56 mV^(-1)

        % Membrane time constant 
        g.taue  = [linspace(4,13,N-1)';20]*1*1e-3;   % 4~13,20 msec
        g.taui  = g.taue*2;                          % 8~26,40 msec

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

    g.state_ue1 =zeros(g.nNodes,1);  
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
%--------------------------------------------------------------------------  
function p=initialize_p()
    
    p.resonance.on = 1; % 1:ON, 0:OFF %plasticity (for sustained resonance)
    p.resonance.memoryLength=1000; %time window for covariance
    p.resonance.WEE_rate=0.05;     %learning rate
    p.resonance.WEI_rate=0.05;     %learning rate
    if p.resonance.on
        disp('Plasticity (bank of oscillators) ON'); 
    else
        disp('Plasticity (bank of oscillators) OFF'); 
    end
end
%--------------------------------------------------------------------------  
function [WEE,WIE,WEI,WII,WEX,WIX,gwE,gwI]=gen_network(N,g)
% input: 
%   N: number of nodes
%  
% output:
%   WEE [NxN]                  ****no negative values***
%   WIE [NxN] to I from E      ****no negative values***
%   WEI [NxN] to E from I      ****no negative values***
%   WII [NxN]                  ****no negative values***
%   gwE [NxN] gaussian weight mask on WEE 
%   gwI [NxN] gaussian weight mask on WEI
    ifirstLayer=[1:N-1]; 
    isecondLayer=N;      
    %----connection strength (node)----
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
    cWE1E2=g.c(5)/length(isecondLayer);
    cWI1E2=g.c(6)/length(isecondLayer);
    cWE1I2=g.c(7)/length(isecondLayer);
    cWI1I2=g.c(8)/length(isecondLayer);

    %------------WEE------------[NxN] to E from E 
    WEE(ifirstLayer,ifirstLayer)=0;
    WEE(isecondLayer,ifirstLayer)=cWE2E1; 
    WEE(ifirstLayer,isecondLayer)=cWE1E2; 
    WEE(eye(N)==1)=cWEE;
    %------------WIE------------[NxN] to I from E 
    WIE(ifirstLayer,ifirstLayer)=0;
    WIE(isecondLayer,ifirstLayer)=cWI2E1; 
    WIE(ifirstLayer,isecondLayer)=cWI1E2; 
    WIE(eye(N)==1)=cWIE;
    %------------WEI------------[NxN] to E from I 
    WEI(ifirstLayer,ifirstLayer)=0;
    WEI(isecondLayer,ifirstLayer)=cWE2I1; 
    WEI(ifirstLayer,isecondLayer)=cWE1I2; 
    WEI(eye(N)==1)=cWEI;
    %------------WII------------[NxN] to I from I 
    WII(ifirstLayer,ifirstLayer)=0;
    WII(isecondLayer,ifirstLayer)=cWI2I1; 
    WII(ifirstLayer,isecondLayer)=cWI1I2; 
    WII(eye(N)==1)=cWII;
    %------------WEX------------[Nx1] to E from input 
    WEX=(220/5)*[ones(N-1,1);0];      
    %------------WIX------------[Nx1] to I from input 
    WIX=WEX*0.5;  
 
    
    %----gaussian weight masks---- 
    gwE=zeros(N,N);
    gwI=zeros(N,N);
    sigmaE=4;
    sigmaI=8;
    for i=ifirstLayer
        x=abs(ifirstLayer-i); %circular
        tmpE=exp(-x.^2/(2*sigmaE^2)); 
        tmpI=exp(-x.^2/(2*sigmaI^2)); 
        gwE(i,ifirstLayer)=tmpE/max(tmpE);
        gwI(i,ifirstLayer)=tmpI/max(tmpI);
    end

end
%--------------------------------------------------------------------------   
% v: potential
% m: firing rate
% v0    = 6e-3;
% e0    = 2.5;
% r     = 560;
function m=sigm(v,v0,e0,r)
  m = (2*e0) ./ (1+exp(r*(v0-v)));
end

