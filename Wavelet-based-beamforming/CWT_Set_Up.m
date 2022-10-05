%  Wavelet based beamforming test code, using SPIRAL array
%  Ver 0.1, working version, 2018 April 9.
%  By Wangqiao Chen & Prof. Xun Huang, LNA, Peking University.
%
% This code can process any mat data file from experiments.
%
% Based on Wavelet & Classic beamforming method
% suit for high speed moving source & unstable source
%
% written by Wangqiao Chen 2018 April 9
%**********************************************************************
%************************** CWT Set Up ********************************
clear all;close all;clc;
C_0 = 340;              % Speed of sound.
sample_fre = 40000;     % Sampling frequency Hz
Num = 62;               % Number of sensors
T_begin = floor(0.3*sample_fre);          % Begining time.
T_during = floor(0.1*sample_fre);         % Beamforming time.
Num_c = 200;            % interpolation number;
main_root = 'D:\MATLAB\CWT_fust\Wavelet_beamforming';
data_path = '\Data';
cwt_path = '\CWT';
data_mat = [main_root,data_path,'\rotating.mat'];
CWT_math = [main_root,cwt_path,'\CWT_rotating.mat'];
% sensor coordinates
mics_coor_mat = ...
    [main_root,'\mics_coor_Array.mat'];
load(mics_coor_mat);          % Array coordinates --> coor
load(data_mat);
longest_r = sqrt(2*(max(coor(:,1))-min(coor(:,1)))^2+...
    2*(max(coor(:,1))-min(coor(:,1)))^2+2*max(coor(:,3))^2);
max_delay = floor(longest_r/C_0*sample_fre);
            % calculate max time delay

T_end = T_begin+T_during+max_delay;     % Stop imaging time
%-------------------------------------------------------------
P_cwt = zeros(Num_c,T_end-T_begin+1,Num);
[test_data,f] = cwt(data(1,T_begin:T_end),sample_fre);
F_cwt = imresize(f,[Num_c,1],'bicubic');  % interpolation for frequency
%------------------------------------------------------------
% Begining CWT
for n = 1:Num
    p_cwt = cwt(data(n,T_begin:T_end),sample_fre);
    for nt=1:T_end-T_begin+1
        P_cwt(:,nt,n) = imresize(p_cwt(:,nt),[Num_c,1],'bicubic');
    end
    clc;disp(['CWT calculated   ',num2str(n),'  /  ',num2str(Num)]);
end
save F_cwt F_cwt
save(CWT_math,'P_cwt','-v7.3');