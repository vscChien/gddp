function [mu,mv,u,v,R]=simulator_3nodes(input,g0)

N=3; % number of nodes
g=initialize_g(N);

if isfield(g0,'ratio1') && isfield(g0,'ratio2') 
    g.ratio1=g0.ratio1; 
    g.ratio2=g0.ratio2; 
end
if isfield(g0,'condition') 
    g.condition=g0.condition; 
end
if isfield(g0,'c31') && isfield(g0,'c13') && isfield(g0,'c32') && isfield(g0,'c23') && isfield(g0,'c21') && isfield(g0,'c12')
    g.c31=g0.c31;
    g.c13=g0.c13;
    g.c32=g0.c32;
    g.c23=g0.c23;         
    g.c21=g0.c21;
    g.c12=g0.c12;
else
    disp('Inter-column connections need to be assigned.')
    return
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
    % ----------plasticity (synaptic adaptation on WEE)----------
    if p.adaptation.on
      da=(1-a)/p.adaptation.taua - p.adaptation.k*a.*repmat(mu(:,t+1)',N,1);
      a=a+da; 
    end 

end 
%----------simulated MEG-----------
MEGweight=([1 1 6]);
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
        
        
        % Default condition
        g.condition=1;
        
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
function p=initialize_p(g)
    % WEE adaptation
    if g.condition==4
      p.adaptation.on = 1;  % 1:ON, 0:OFF 
    else
      p.adaptation.on = 0;     
    end
    p.adaptation.taua=200;  % 200 msec
    p.adaptation.k=2*1e-3;  % decay rate

    if p.adaptation.on, disp('WEE adaptation ON'); else, disp('WEE adaptation OFF'); end
end
%--------------------------------------------------------------------------  
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
%----connection strength (node)----
cWEE=0.8;
cWIE=0.6;
cWEI=0.2;           
cWII=0.05;    

%------------WEE------------[NxN] to E from E 
WEE=[cWEE     g.c12(1) g.c13(1);...  
     g.c21(1) cWEE     g.c23(1);...
     g.c31(1) g.c32(1) cWEE   ];
a=ones(N); % adaptation term 
%------------WEI------------[NxN] to I from E 
WIE=[cWIE     g.c12(2) g.c13(2);...  
     g.c21(2) cWIE     g.c23(2);... 
     g.c31(2) g.c32(2) cWIE   ];
%------------WIE------------[NxN] to E from I 
WEI=[cWEI     g.c12(3) g.c13(3);...  
     g.c21(3) cWEI     g.c23(3);... 
     g.c31(3) g.c32(3) cWEI   ];
%------------WII------------[NxN] to I from I 
WII=[cWII     g.c12(4) g.c13(4);...  
     g.c21(4) cWII     g.c23(4);... 
     g.c31(4) g.c32(4) cWII   ];
%------------WEX------------[Nx2] to E from inputs 
WEX=(220/5)*[1,0;...
             0,1;...
             0,0];      
%------------WIX------------[Nx2] to I from inputs 
WIX=WEX*0.5; 

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