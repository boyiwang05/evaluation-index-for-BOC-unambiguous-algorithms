
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
BW=30*1.023e6;d=0.1*Tc;%BOC(10,5)相关器间隔
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_BW=10000;
f=linspace(-BW/2,BW/2,N_BW);
PSD_BOC=PSDcal_BOCs(f, fs, Tc);
power_loss_filter_dB=10*log10(trapz(f,PSD_BOC));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=1;%载波功率

N_C=1;
% C_N0_dB=linspace(20,45,N_C)-power_loss_filter_dB;
C_N0_dB=25;
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
n_loop=100000;

Trackerror1=zeros(1,n_loop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%BOC(14,2)%%%%%%%%%
% Dz=1.652e-4;
Dz=5.27e-4;
%%
error_x1=zeros(1,n_loop);
error_y1=zeros(1,n_loop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MDR_dB=-6;
MDR=10^(MDR_dB/20);
% N_multipath=40;
N_multipath=1;
delay_t=0.4*Tc;
% delay_t=linspace(0,1.2*Tc,N_multipath);
delay_m=delay_t*3e8;
% Multipath_error=zeros(2,N_multipath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S_CA_receiver,t]=yt_BOCs_function(c,t_begin,t_end,Tc,fs,f_sample);
N_zong=length(t);

TrackErrRMS=zeros(1,N_C);
h_wait=waitbar(0);

% Num_false_acquisition=200;
% t_ref_begin_scan=linspace(-1*Tc,Tc,Num_false_acquisition);

for m=1:N_multipath
    for k=1:N_C
        tic;
        t0=0.0*Ts;%真实传播时间，初始误差在一个码片内

        t_begin=t0;
        t_end=t_begin+Tp;
%         t_ref_begin=t_ref_begin_scan(m);
        t_ref_begin1=t0;
        t_ref_begin2=t0;
        
        for n=1:n_loop

            t_ref_end1=t_ref_begin1+Tp;t_ref_end2=t_ref_begin2+Tp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% 0 phase loop signal generator
         [s_ref1_E_1,s_ref2_E_1,t_E_1]=yt_WPC_reference_waveform_function(c,t_ref_begin1+d/2,t_ref_end1+d/2,Tc,fs,f_sample);
         [s_ref1_L_1,s_ref2_L_1,t_L_1]=yt_WPC_reference_waveform_function(c,t_ref_begin1-d/2,t_ref_end1-d/2,Tc,fs,f_sample);
         [s_ref1_P_1,s_ref2_P_1,t_P_1]=yt_WPC_reference_waveform_function(c,t_ref_begin1,t_ref_end1,Tc,fs,f_sample);    
         %%   180 phase loop signal generator
         [s_ref1_E_2,s_ref2_E_2,t_E_2]=yt_WPC_reference_waveform_function(c,t_ref_begin2+d/2,t_ref_end2+d/2,Tc,fs,f_sample);
         [s_ref1_L_2,s_ref2_L_2,t_L_2]=yt_WPC_reference_waveform_function(c,t_ref_begin2-d/2,t_ref_end2-d/2,Tc,fs,f_sample);
         [s_ref1_P_2,s_ref2_P_2,t_P_2]=yt_WPC_reference_waveform_function(c,t_ref_begin2,t_ref_end2,Tc,fs,f_sample);   

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             N_zong=length(t);
             %%%%%%%%%%%%%%%%%%%%%%%%%%   
             %%
            S_receiver_I=sqrt(2*C)*S_CA_receiver+sqrt(I_NoisePower(k))*randn(1,N_zong);
            S_receiver_Q=sqrt(Q_NoisePower(k))*randn(1,N_zong);
                     %%%%%%%%%%%%%%%%%%%%%%%%%%
         %多径信号

         [S_CA_receiver_multipath,t]=yt_BOCs_function(c,t_begin-delay_t(m),t_end-delay_t(m),Tc,fs,f_sample);
         S_receiver_1=(S_receiver_I+1i*S_receiver_Q)+MDR*S_CA_receiver_multipath;        
         S_receiver_2=(S_receiver_I+1i*S_receiver_Q)-MDR*S_CA_receiver_multipath;

            fft_sig  = fft(S_receiver_1);
            L_res   = round(BW/f_sample/2*N_zong);
            fft_sig(L_res+1:end-L_res)  = 0;
            sigBandL = ifft(fft_sig);
            RecSigI1 = real(sigBandL);
            RecSigQ1 = imag(sigBandL);
            
            fft_sig  = fft(S_receiver_2);
            L_res   = round(BW/f_sample/2*N_zong);
            fft_sig(L_res+1:end-L_res)  = 0;
            sigBandL = ifft(fft_sig);
            RecSigI2 = real(sigBandL);
            RecSigQ2 = imag(sigBandL);

            
            %% 0 phase multipath
            IE1 = sum(s_ref1_E_1 .* RecSigI1)/N_zong;IE2 = sum(s_ref2_E_1 .* RecSigI1)/N_zong;
            IL1 = sum(s_ref1_L_1.* RecSigI1)/N_zong;IL2 = sum(s_ref2_L_1.* RecSigI1)/N_zong;
            QE1 = sum(s_ref1_E_1.* RecSigQ1)/N_zong;QE2 = sum(s_ref2_E_1.* RecSigQ1)/N_zong;
            QL1 = sum(s_ref1_L_1.* RecSigQ1)/N_zong;QL2= sum(s_ref2_L_1.* RecSigQ1)/N_zong;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %% 鉴别器1
            DiscrimOut1  =IE1*IE2+QE1*QE2-IL1*IL2-QL1*QL2;%鉴别器1
            
             %% 180 phase multipath
            IE1 = sum(s_ref1_E_2 .* RecSigI2)/N_zong;IE2 = sum(s_ref2_E_2 .* RecSigI2)/N_zong;
            IL1 = sum(s_ref1_L_2.* RecSigI2)/N_zong;IL2 = sum(s_ref2_L_2.* RecSigI2)/N_zong;
            QE1 = sum(s_ref1_E_2.* RecSigQ2)/N_zong;QE2 = sum(s_ref2_E_2.* RecSigQ2)/N_zong;
            QL1 = sum(s_ref1_L_2.* RecSigQ2)/N_zong;QL2= sum(s_ref2_L_2.* RecSigQ2)/N_zong;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %% 鉴别器1
            DiscrimOut2  =IE1*IE2+QE1*QE2-IL1*IL2-QL1*QL2;%鉴别器1

        %% 0 phase loop update
            error_x1(n)=(t_begin-t_ref_begin1)/Tc;%真实误差
            error_y1(n)=DiscrimOut1;%鉴别器输出


            Filter_output1=DiscrimOut1*Dz;%一阶环路滤波器

            t_ref_begin1=t_ref_begin1+Filter_output1*Tc;

            Trackerror1(n)=(t_ref_begin1-t_begin)*3e8;%跟踪误差
            
            %% 180 phase loop update
            error_x2(n)=(t_begin-t_ref_begin2)/Tc;%真实误差
            error_y2(n)=DiscrimOut2;%鉴别器输出


            Filter_output2=DiscrimOut2*Dz;%一阶环路滤波器

            t_ref_begin2=t_ref_begin2+Filter_output2*Tc;

            Trackerror1(n)=(t_ref_begin2-t_begin)*3e8;%跟踪误差
            %%

            temp_string=['已运行' num2str(floor(((k-1)*n_loop+n+m)/(N_C*n_loop*N_multipath)*10000)/100) '%'];
            waitbar(((k-1)*n_loop+n+m)/(N_C*n_loop*N_multipath),h_wait,temp_string);
        end
    end
    TrackErrSTD_0_phase(m)=std(error_x1);
    TrackErrSTD_180_phase(m)=std(error_x2);
    
    TrackErrMEAN_0_phase(m)=mean(error_x1);
    TrackErrMEAN_180_phase(m)=mean(error_x2);
    
    TrackErr_0_phase(m,:)=error_x1;
    TrackErr_180_phase(m,:)=error_x2;
    
    %% false lock counter
%     false_lock_counter=0;
%     if mean(error_x1)>0.9*Ts || mean(error_x2)>0.9*Ts
%         false_lock_counter=false_lock_counter+1;
%         
%     end
    
    toc;
end
savefile='multipath_induced_bias_yans_method_BOC_10_5_30M_01Tc_25dBHz.mat';
save(savefile);


close(h_wait);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot false lock index
% figure;plot(delay_t/Tc,TrackErrSTD_0_phase,'Linewidth',2);grid on;
% xlabel('multipath delay (chips)');
% ylabel('code tracking error std (chips)');
% saveas(gcf,'0_phase_multipath_induced_bias_yans_method_BOC_10_5_30M_01Tc_45dBHz');
% 
% figure;plot(delay_t/Tc,TrackErrSTD_180_phase,'Linewidth',2);grid on;
% xlabel('multipath delay (chips)');
% ylabel('code tracking error std (chips)');
% saveas(gcf,'180_phase_multipath_induced_bias_yans_method_BOC_10_5_30M_01Tc_45dBHz');