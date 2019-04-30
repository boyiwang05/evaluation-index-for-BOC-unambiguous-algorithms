function [s1,s2,t]=yt_BOCs_SPAR_function(c,t_begin,t_end,Tc,fs,fsample,Pulse_width)
%c是一个码周期的扩频码序列
%t_begin起始时间
%t_end终止时间
%Tc码片宽度
%fs子载波频率
%fsample采样频率
% M=2*fs*Tc;
L_code=length(c);
% L_SUB=M;
% for m=1:(M/2-1)
%     d(m)=(-1)^(m-1);
%     d(M-1-m)=(-1)^(m-1);
% end
% d(2:M-1)=d(1:M-2);
% d(1)=M-1;
% d(M)=M-1;
% d=d/sqrt(2*M-1);
t=t_begin:1/fsample:t_end;
N_Tc=mod(floor(t/Tc),L_code)+1;%每个时间点对应的码
code=c(N_Tc);
N_Sub=mod(t,Tc);
for i=1:length(t)
    if N_Sub(i)<Pulse_width*Tc
        s1(i)=1;
        s2(i)=1;
    elseif N_Sub(i)>(1-Pulse_width)*Tc
        s1(i)=-1;
        s2(i)=1;   
    else 
        s1(i)=0;
        s2(i)=0;
    end
end
s1=code.*s1/2/Pulse_width;
s2=code.*s2/2/Pulse_width;