function [inputE,inputI,times]=gen_input (N,SOA,onofftime,input1ratio,input2ratio)

if input1ratio>1 | input1ratio<0
    error('0<=input1ratio<=1')
end

ifirstLayer=1; % index of first-layer nodes
isecondLayer=2; % index of second-layer nodes 

dt              = 1e-3;   % 1000 Hz
times           = dt:dt:7; % 1 sec 

background=2.5; 
amplitude=1;  
JIratio=0.5;


input=zeros(N,length(times));% N x times
         if SOA==0 % --- square ---       
              offset=1000;%ms
              code.type=4; % pure-tone sequence (method: piecewise)
              code.n=1;
              code.tlength=length(times);%1000ms
              code.dt=0.001; %1ms
              code.peak=4-background;%pps
              code.soa=10000;%ms
              code.duration=diff(onofftime);%3000;%ms
              code.risetime=10;
              code.offset=onofftime(1);
              code.figureon=0;

              x1=gen_input_core(code)'; %x1(end-offset:end)=0; 
              %x1(3900:4100)=0; %omission
              x=[x1];
              inputE= repmat(background+x*amplitude*input1ratio        ,N,1); 
              inputI= repmat(           x*amplitude*input1ratio*JIratio,N,1); 
         elseif SOA~=0 %-----square periodic----
              %offset=2050;%ms
              JIratio=0.5;%0.1;%0.05;  
              code.type=4; % pure-tone sequence (method: piecewise)
              code.n=1;
              code.tlength=length(times);%1000ms
              code.dt=0.001; %1ms
              code.peak=4-background;%pps
              code.soa=SOA;%400;%ms
              code.duration=100;%50;%ms
              code.risetime=10;
              code.offset=onofftime(1);
              code.figureon=0;

              nCycles=floor(diff(onofftime)/code.soa);
              x1=gen_input_core(code)'; x1(4000+code.soa*nCycles:end)=0; 
 
              x=[x1];
              inputE= repmat(background+x*amplitude*input1ratio        ,N,1); 
              inputI= repmat(           x*amplitude*input1ratio*JIratio,N,1); 
         end
%          xx=1.2*(1-input1ratio-0.1)/(1-0.1);
%          xx(xx<0)=0;
         inputE(isecondLayer,:)=background+x*amplitude*input2ratio; % input to second layer
         inputI(isecondLayer,:)=x*amplitude*JIratio*input2ratio; % input to second layer
         
%          inputE(isecondLayer,:)=background; % input to second layer
%          inputI(isecondLayer,:)=0; % input to second layer
              
end
