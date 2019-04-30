
clear;
clc;
%%
%正弦BOC参数
c=CA_code(1);%得到CA码序列1
L_CA=length(c);%CA码序列长度
m=10;n=5;
Rc=n*1.023e6;%码速率
Tc=1/Rc;%码片长度
f_sample=100e6;%采样频率
T_sample=1/f_sample;
Tp=1e-3-T_sample;%相干积分时间1ms
fs=m*1.023e6;
Ts=1/fs/2;
%%
%%%相关器间隔
BW=30*1.023e6;d_Subcarrier=0.1*Tc;d_code=Ts;%BOC(10,5)相关器间隔%shiyan2
Dz_code  = 7.509e-4;%K_code=5.3129
Dz_Subcarrier  =1.3331e-4;%K_Subcarrier=29.9243
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_BW=10000;
f=linspace(-BW/2,BW/2,N_BW);
PSD_BOC=PSDcal_BOCs(f, fs, Tc);
power_loss_filter_dB=10*log10(trapz(f,PSD_BOC));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=1;%载波功率

N_C=1;
% C_N0_dB=linspace(20,45,N_C)-power_loss_filter_dB;
C_N0_dB=45;
C_N0=10.^(C_N0_dB/10);
N0=C./C_N0;%单边噪声功率谱密度
I_NoisePower=N0*BW*f_sample/BW;%在接收带宽内的噪声功率,因为在滤波后，噪声功率减为此处的BW/fsample，为了使滤波后噪声功率保持2N0*BW(射频变基带噪声谱密度变2倍)，在此处乘以fsample/BW
Q_NoisePower=N0*f_sample;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=0.0*Ts;%真实传播时间，初始误差在一个码片内
% t0=-0.9648*Tc;%真introduce false lock
t_begin=t0;
t_end=t_begin+Tp;
t_ref_begin_scan=0;
% n_loop=100000;
n_loop=10000;

real_error=zeros(1,n_loop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
error_x=zeros(1,n_loop);
error_y=zeros(1,n_loop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S_CA_receiver,t]=yt_BOCs_function(c,t_begin,t_end,Tc,fs,f_sample);
N_zong=length(t);

TrackErrRMS=zeros(1,N_C);
h_wait=waitbar(0);

Num_false_acquisition=200;
t_ref_begin_scan=linspace(-1*Tc,Tc,Num_false_acquisition);

for m=1:Num_false_acquisition
    for k=1:N_C
        tic;
        t0=0.0*Ts;%真实传播时间，初始误差在一个码片内

        t_begin=t0;
        t_end=t_begin+Tp;
        t_code_ref_begin=t_ref_begin_scan(m);
        t_Subcarrier_ref_begin=t_ref_begin_scan(m);
        t_estimate_begin=0;
        
        for n=1:n_loop

        t_code_ref_end=t_code_ref_begin+Tp;t_Subcarrier_ref_end=t_Subcarrier_ref_begin+Tp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    [s_code_ref_E,t]=yt_BOCs_function_Dual_Estimate_Solution(c,t_code_ref_begin+d_code/2,t_code_ref_end+d_code/2,t_Subcarrier_ref_begin,t_Subcarrier_ref_end,Tc,fs,f_sample);
    [s_code_ref_L,t]=yt_BOCs_function_Dual_Estimate_Solution(c,t_code_ref_begin-d_code/2,t_code_ref_end-d_code/2,t_Subcarrier_ref_begin,t_Subcarrier_ref_end,Tc,fs,f_sample);    
    [s_Subcarrier_ref_E,t]=yt_BOCs_function_Dual_Estimate_Solution(c,t_code_ref_begin,t_code_ref_end,t_Subcarrier_ref_begin+d_Subcarrier/2,t_Subcarrier_ref_end+d_Subcarrier/2,Tc,fs,f_sample);
    [s_Subcarrier_ref_L,t]=yt_BOCs_function_Dual_Estimate_Solution(c,t_code_ref_begin,t_code_ref_end,t_Subcarrier_ref_begin-d_Subcarrier/2,t_Subcarrier_ref_end-d_Subcarrier/2,Tc,fs,f_sample);       
     N_zong=length(t);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    滤波
    fft_sig  = fft(sqrt(2*C)*S_CA_receiver+sqrt(I_NoisePower(k))*randn(1,N_zong)+1i*sqrt(Q_NoisePower(k))*randn(1,N_zong));
    L_res   = round(BW/f_sample/2*N_zong);
    fft_sig(L_res+1:end-L_res)  = 0;
    sigBandL = ifft(fft_sig);
    %%%%%%%%%%%%%%   
    RecSigI = real(sigBandL);
    RecSigQ = imag(sigBandL);
    %%
        Int_code_Ie = sum(s_code_ref_E.* RecSigI)/N_zong;Int_code_Il = sum(s_code_ref_L.* RecSigI)/N_zong;
        Int_code_Qe = sum(s_code_ref_E.* RecSigQ)/N_zong;Int_code_Ql = sum(s_code_ref_L.* RecSigQ)/N_zong;
        
        Int_Subcarrier_Ie = sum(s_Subcarrier_ref_E.* RecSigI)/N_zong;Int_Subcarrier_Il = sum(s_Subcarrier_ref_L.* RecSigI)/N_zong;
        Int_Subcarrier_Qe = sum(s_Subcarrier_ref_E.* RecSigQ)/N_zong;Int_Subcarrier_Ql = sum(s_Subcarrier_ref_L.* RecSigQ)/N_zong;
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%非相干鉴别函数    
%%
        DiscrimOut_code= (Int_code_Ie^2+Int_code_Qe^2)-(Int_code_Il^2+Int_code_Ql^2);%
        DiscrimOut_Subcarrier= (Int_Subcarrier_Ie^2+Int_Subcarrier_Qe^2)-(Int_Subcarrier_Il^2+Int_Subcarrier_Ql^2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        error_code_x(n)=(t_begin-t_code_ref_begin)/Tc;%真实误差
        error_code_y(n)=DiscrimOut_code;%鉴别器输出
        error_Subcarrier_x(n)=(t_begin-t_Subcarrier_ref_begin)/Tc;%真实误差
        error_Subcarrier_y(n)=DiscrimOut_Subcarrier;%鉴别器输出
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
        Filter_code_output=DiscrimOut_code*Dz_code;%一阶环路滤波器
        Filter_Subcarrier_output=DiscrimOut_Subcarrier*Dz_Subcarrier;%一阶环路滤波器

       %%
        %更新时延估计
        t_code_ref_begin=t_code_ref_begin+Filter_code_output*Tc;
        t_Subcarrier_ref_begin=t_Subcarrier_ref_begin+Filter_Subcarrier_output*Tc;
       %%
        t_estimate_begin=t_Subcarrier_ref_begin+round((t_code_ref_begin-t_Subcarrier_ref_begin)/Ts)*Ts;%真实估计时间
        real_error(n)=(t_begin-t_estimate_begin)/Tc;%跟踪误差
        
        temp_state(n)=abs(round((t_code_ref_begin-t_Subcarrier_ref_begin)/Ts));

            temp_string=['已运行' num2str(floor(((k-1)*n_loop+n+m)/(N_C*n_loop*Num_false_acquisition)*10000)/100) '%'];
            waitbar(((k-1)*n_loop+n+m)/(N_C*n_loop*Num_false_acquisition),h_wait,temp_string);
        end
    end
    code_TrackErrSTD(m)=std(error_code_x);
    subcarrier_TrackErrSTD(m)=std(error_Subcarrier_x);
    real_TrackErrSTD(m)=std(real_error);
    
    error_code_recorder(m,:)=error_code_x;
    error_Subcarrier_recorder(m,:)=error_Subcarrier_x;
    real_error_recorder(m,:)=real_error;
    
    toc;
end
savefile='false_lock_index_DET_BOC_10_5_30M_01Tc_45dBHz.mat';
save(savefile);


close(h_wait);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot false lock index
figure;plot(t_ref_begin_scan/Tc,code_TrackErrSTD,t_ref_begin_scan/Tc,subcarrier_TrackErrSTD,t_ref_begin_scan/Tc,real_TrackErrSTD,'Linewidth',2);grid on;
xlabel('false acuiqisition results (chips)');
ylabel('code tracking error std (chips)');
legend('code loop','subcarrier loop','real DET');
saveas(gcf,'false_lock_index_DET_BOC_10_5_30M_01Tc_45dBHz');
