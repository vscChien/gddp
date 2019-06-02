% Description:
%   To generate different types of inputs to the network
%   
%
% Inputs:
%     code         : [struct]
%     code.type    : type of input to generate (e.g., background noise, periodic sequence, etc.)
%     code.n       : number of nodes
%     code.tlength : length of input (time points)
%                    e.g.
%                       1000 %1 sec, if code.dt==10^-3
%     code.dt      : sampling interval (sec)
%                    e.g.
%                       1e-3 % 0.001 sec                   
%     code.*       : other parameters for different types 
%
% Outputs:
%     data         : envolope, [code.tlength x code.n] array 
%     timepoints   : [code.tlength x 1] array
%     dataCode     : in order to put different pitch in the envolope [code.tlength x code.n] array 
%
% Example:
%     >> code.type=1;     % generate background noise 
%     >> code.n=5;
%     >> code.tlength=1000;
%     >> code.dt=1e-3;
%     >> code.mu=5;       % mean value 5*44 spikes/sec
%     >> code.sigma=0.1;  % standard deviation 
%     >> code.offset=[100,175,175,175,175];  % the time in ms for each node to kick off.
%     >> code.figureon=1; % show figure
%     >> [data,timepoints]=gen_input_core(code);
%
% 20160713 Vincent
function [data,timepoints,dataCode]=gen_input_core(code)
data=[];
dataCode=[];
timepoints=[];

if ~isstruct(code)
    disp('-----------------')
    disp(mfilename)
    disp('Input ''code'' should be a structure: Please see description.')
    disp('-----------------')
    return
end

%----generate data----
switch code.type
    case 1 % generate background noise
        data=gen_input_noise(code.n,...
                             code.tlength,...
                             code.mu,...
                             code.sigma,...
                             code.offset);
    case 2 % generate oscillating sequence (method: kernal)
        data=gen_input_osci(code.n,...
                            code.tlength,...
                            code.peak,...
                            code.dt,...
                            code.freq,...
                            code.offset,...
                            code.reverseSequence);
    case 3 % generate oscillating sequence (method: piecewise)
        data=gen_input_osci2(code.n,...
                            code.tlength,...
                            code.peak,...
                            code.dt,...
                            code.freq,...
                            code.risetime,...
                            code.offset,...
                            code.reverseSequence);
    case 4 %generate pure-tone sequence (method: piecewise)
       code.n=1;
       data=gen_input_tone(code.tlength,...
                            code.dt,...
                            code.peak,...
                            code.soa,...
                            code.duration,...% tone ON duration
                            code.risetime,...
                            code.offset);  % start time       
    case 5 %generate two-tone pair sequence (method: piecewise)
           % Note: 1st tone(+), 2nd tone(-)
       code.n=1;
       data=gen_input_tonePair(code.tlength,...
                            code.dt,...
                            code.peak,...
                            code.soa,...
                            code.isi,...
                            code.duration,...% tone ON duration
                            code.risetime,...
                            code.offset);  % start time   
    case 6 %generate N-tone sequence in sequence (method: piecewise)
       code.n=1;
       [data,dataCode]=gen_input_NtoneSeq(code.tlength,...
                               code.dt,...
                               code.peak,...
                               code.seqDuration,...%code.soa,...
                               code.seqSOA,...
                               code.freqCode,...
                               code.duration,...% tone ON duration
                               code.risetime,...
                               code.offset);  % start time  
    case 7 %generate SD mask with p(S)=0.9 in a sequence (method: piecewise)
       code.n=1;
       [data]=gen_input_SDmask(code.tlength,...
                               code.dt,...
                               code.soa,...
                               code.offset,...
                               code.p);  % start time                        
end
timepoints=[1:size(data,1)]'*code.dt;

%----show figure----
if code.figureon==1
  figure;plot(timepoints,data);
  xlabel('sec')
  ylabel('mean firing rate (*44 spikes/sec)')
  legend(num2str([1:code.n]'))
end


%==========================================================================
% generate background noise
% Example:
%     >> code.type=1;        % background noise 
%     >> code.n=5;           % number of nodes
%     >> code.tlength=1000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.mu=5;          % mean value of the background noise (5*44 spikes/sec)
%     >> code.sigma=0.1;     % standard deviation 
%     >> code.offset=[100,175,175,175,175];  % time points for each node to kick off (being not 0).
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);
function data=gen_input_noise(n,tlength,mu,sigma,offset)
    data=ones(tlength,n)*mu + randn(tlength,n)*sigma;
    for i=1:n
      data(1:1+offset(i),i)=0; 
    end 
    data(data<0)=0;
    
    
    

%==========================================================================
% generate oscillating sequence
% Example:
%     >> code.type=2;        % oscillating sequence
%     >> code.n=5;           % number of nodes
%     >> code.tlength=1000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.peak=5;        % peak value of the sequence (5->5*44 spikes/sec)
%     >> code.freq=8;        % periodicity of the sequence (Hz), 
%                               ex. freq=8 means 8 peaks (each input) in a sec  
%     >> code.offset=100;    % start the sequence at 100 time points.
%     >> code.reverseSequence=0; %from 1-2-3-4 to 4-3-2-1
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);
function data=gen_input_osci(n,tlength,peak,dt,freq,offset,reverseSequence)
    %------------kernel--------------
    switch 2
        case 1 % kernel 1: normal distribution     
            x=-100:1:100; % ms
            mu=0; sig=1; % mean; std
            y=normpdf(x,mu,sig);  %y=exp(-(x-mu).^2./(2*sig^2))/(sig*sqrt(2*pi));
            y=y*peak/max(y); % normalize to 5
        case 2  % kernel 2: bell, http://de.mathworks.com/help/fuzzy/gbellmf.html    
            x=-50:1:50; % ms
            a=15; % width
            b=5; % steep
            c=0; % mean
            y=1./(1+abs((x-c)./a).^(2*b));
            y=y*peak/max(y); % normalize to 5
    end
    %-----------------
    data=[];
    for i=1:n
      tt=zeros(1,tlength);
      tt(round(1+offset+(i-1)*(1/dt/freq/n):1/dt/freq:end))=1;
      data(i,:)=conv(tt,y,'same');
    end
    data=data';
    %-------reverse sequence or not------
    if reverseSequence  %from 1-2-3-4 to 4-3-2-1
      data=fliplr(data);
    end

    
%==========================================================================
% generate oscillating sequence
% Example:
%     >> code.type=3;        % oscillating sequence
%     >> code.n=5;           % number of nodes
%     >> code.tlength=4000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.peak=5;        % peak value of the sequence (5->5*44 spikes/sec)
%     >> code.risetime=5;    % rise time of the kernel (should be less than a period)
%     >> code.freq=8;        % periodicity of the sequence (Hz), 
%                               ex. freq=8 means 8 peaks (each input) in a sec  
%     >> code.offset=100;    % start the sequence at 100 time points.
%     >> code.reverseSequence=0; %from 1-2-3-4 to 4-3-2-1
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);
function data=gen_input_osci2(n,tlength,peak,dt,freq,risetime,offset,reverseSequence)
    %------------kernel--------------
    inputSeqPeriod=round(1/(freq*n*dt)); %250ms
    if inputSeqPeriod<=risetime
        error('code.risetime should be less than 1/(code.freq)')
    end
    duration=inputSeqPeriod+risetime;%ms 
    u1=interp1([1,risetime+1,duration-risetime,duration],[0,1,1,0]*peak,1:duration); % kernel: 50 ms duration (5ms rise and fall)
    %-----------------
    data=zeros(tlength,n);
    for i=1:n
        startpoints=[1+(i-1)*inputSeqPeriod+offset:inputSeqPeriod*n:tlength];
        for ss=startpoints
            if ss+length(u1)-1<=length(data)
               data(ss:ss+length(u1)-1,i)=u1;
            end
        end
    end
    %-------reverse sequence or not------
    if reverseSequence  %from 1-2-3-4 to 4-3-2-1
      data=fliplr(data);
    end
    
    
%==========================================================================
% generate pure-tone sequence
% Example:
%     >> code.type=4;        % pure-tone sequence
%     >> code.tlength=4000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.peak=5;        % peak value of the sequence (5->5*44 spikes/sec)
%     >> code.soa=200;       % msec
%     >> code.duration=50;   % tone ON duration (msec) including rising/falling time
%     >> code.risetime=5;    % rise time of the kernel (should be less than a period)
%     >> code.offset=100;    % start the sequence at 100 time points.
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);
    
function   data=gen_input_tone(tlength,dt,peak,soa,duration,risetime,offset)  
    % msec to timepoints
    tlength=tlength/dt/1000;
    soa=soa/dt/1000;
    duration=duration/dt/1000;
    risetime=risetime/dt/1000;
    offset=round(offset/dt/1000);
    %
    u1=interp1([0,risetime,duration-risetime,duration,soa],[0,1,1,0,0]*peak,1:soa); % kernel: 50 ms duration (5ms rise and fall)
    data=repmat(u1,1,ceil(tlength/soa));                        
    data=[zeros(1,offset),data];
    data=data(1:tlength);
    data=data';
    
%==========================================================================
% generate pure-tone pair sequence
% Example:
%     >> code.type=5;        % pure-tone sequence
%     >> code.tlength=4000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.peak=5;        % peak value of the sequence (5->5*44 spikes/sec)
%     >> code.soa=500;       % msec SOA between pairs
%     >> code.ISI=100;       % msec ISI in the pair
%     >> code.duration=50;   % tone ON duration (msec) including rising/falling time
%     >> code.risetime=5;    % rise time of the kernel (should be less than a period)
%     >> code.offset=100;    % start the sequence at 100 time points.
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);
    
function   data=gen_input_tonePair(tlength,dt,peak,soa,isi,duration,risetime,offset)  
    % msec to timepoints
    tlength=tlength/dt/1000;
    soa=soa/dt/1000;
    isi=isi/dt/1000;
    duration=duration/dt/1000;
    risetime=risetime/dt/1000;
    offset=round(offset/dt/1000);
    %
    u1=interp1([1,risetime+1,duration-risetime,duration,...
                duration+isi,duration+isi+risetime,duration*2+isi-risetime,duration*2+isi...
                soa],[0,1,1,0,0,-1,-1,0,0]*peak,1:soa); % kernel: 50 ms duration (5ms rise and fall)
    data=repmat(u1,1,ceil(tlength/soa));                        
    data=[zeros(1,offset),data];
    data=data(1:tlength);
    data=data'; 
    
    
    
%==========================================================================
% generate N-tone sequence
% Example:
%     >> code.type=6;        % N-tone sequence in a sequence
%     >> code.tlength=4000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.peak=5;        % peak value of the sequence (5->5*44 spikes/sec)
%     >> code.seqSOA=500*5;
%     >> code.freqCode=[1000,1100,1200,1300,1400]; % N-tone sequence 
%     >> code.duration=50;   % tone ON duration (msec) including rising/falling time
%     >> code.risetime=5;    % rise time of the kernel (should be less than a period)
%     >> code.offset=100;    % start the sequence at 100 time points.
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);    
%function   [data,dataCode]=gen_input_NtoneSeq(tlength,dt,peak,seqDuration,soa,seqSOA,freqCode,duration,risetime,offset)  
function   [data,dataCode]=gen_input_NtoneSeq(tlength,dt,peak,seqDuration,seqSOA,freqCode,duration,risetime,offset)  

%     if seqDuration > seqSOA
%        error('sequence length should be < SOA!') 
%     end

    nTones=length(freqCode);
    if mod(seqSOA,nTones)~=0
        error('SOA should be multiple of number of tones in a sequence!')
    else
        toneSOA=seqSOA/nTones;
    end
    tlength=tlength/dt/1000;
    toneSOA=toneSOA/dt/1000;
    duration=duration/dt/1000;
    risetime=risetime/dt/1000;
    offset=round(offset/dt/1000);
    %
    %u1=interp1([1,risetime+1,duration-risetime,duration,toneSOA],[0,1,1,0,0]*peak,1:toneSOA); % kernel: 50 ms duration (5ms rise and fall)
    u1=interp1([1,risetime+1,duration-risetime,duration],[0,1,1,0]*peak,1:duration); % kernel: 50 ms duration (5ms rise and fall)
    u=(u1>0);
    u1=repmat(u1,1,nTones);
    u1Code=[];
    for i=1:nTones
        u1Code=[u1Code,u*freqCode(i)];
    end
    u1=[u1,zeros(1,seqSOA-length(u1))];
    u1Code=[u1Code,zeros(1,seqSOA-length(u1Code))];
    
    data=repmat(u1,1,ceil(tlength/seqSOA));  
    dataCode=repmat(u1Code,1,ceil(tlength/seqSOA));
    data=[zeros(1,offset),data];
    dataCode=[zeros(1,offset),dataCode];
    data=data(1:tlength);
    dataCode=dataCode(1:tlength);
    data=data'; 
    dataCode=dataCode';
        
%==========================================================================
% generate SD mask with p(S)=0.9 in a sequence (method: piecewise)
% Example:
%     >> code.type=7;        % N-tone sequence in a sequence
%     >> code.tlength=4000;  % time points
%     >> code.dt=1e-3;       % sampling interval (sec)
%     >> code.soa=500;       % msec SOA
%     >> code.offset=100;    % start the sequence at 100 time points.
%     >> code.p=0.9;         % probability of S (standard)
%     >> code.figureon=1;    % show figure
%     >> [data,timepoints]=gen_input_core(code);          
function   [data]=gen_input_SDmask(tlength,dt,soa,offset,p)     
    % msec to timepoints
    tlength=tlength/dt/1000;
    soa=soa/dt/1000;  
    offset=round(offset/dt/1000);
    %----------------
    nRepetitions=ceil(tlength/soa/p);
    mask=[]; % [nRepetitions x p]
    for i=1:nRepetitions
        tmp=[0,ones(1,p-2)]; %[1x(p-1)]
        tmp=[1,tmp(randperm(length(tmp)))]; %[1xp] where 1st will always be 1
        mask=[mask;tmp]; 
    end
    mask=mask'; % [p x nRepetitions]
    mask=mask(:)'; % [1x (p*nRepetitions)]
    mask=repmat(mask,soa,1); % [soa x (p*nRepetitions)]
    mask=mask(:)';
    %-----------                      
    data=[zeros(1,offset),mask];
    data=data(1:tlength);
    data=data'; 
        
        
        
