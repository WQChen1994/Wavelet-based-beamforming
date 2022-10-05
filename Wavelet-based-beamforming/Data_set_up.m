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
%************************** Data Set Up********************************
clear all;close all;clc;
C_0 = 340;              % Speed of sound.
sample_fre = 40000;     % Sampling frequency Hz
over_sample = 25;       % over sample radio
over_sample_fre = over_sample*sample_fre; % over sample frequency Hz
Fre = 2000;             % signal frequency
Num = 62;               % Number of sensors
t=linspace(0,1,over_sample_fre);
%-----------------------------------------------------
% root
main_root = 'D:\MATLAB\CWT_fust\Wavelet_beamforming';
data_path = '\Data';
data_mat = [main_root,data_path,'\rotating.mat'];
% sensor coordinates
mics_coor_mat = ...
    [main_root,'\mics_coor_Array.mat'];
load(mics_coor_mat);          % Array coordinates --> coor
%----------------------------------------------------
%************************* Moving motion ******************************
% Rotation center at 0,0
Rot_fre = 50;           % Hz;
Rot_R  = 0.3;           % m;
Rot_v = Rot_R*Rot_fre*2*pi; % rotating speed.
phi = Rot_fre*2*pi*t;   % Rotating phase.
Rot_x = Rot_R*cos(phi);
Rot_y = -Rot_R*sin(phi); % Source location.
Rot_vx = Rot_y.*Rot_fre*2*pi;
Rot_vy = -Rot_x.*Rot_fre*2*pi; % Source speed.
%------------------------------------------------------
% Signal generate
d_x = zeros(over_sample_fre,1);
d_y = zeros(over_sample_fre,1);
r = zeros(over_sample_fre,1);
alpha = zeros(over_sample_fre,1);
beta = zeros(over_sample_fre,1);
d_t = zeros(Num,over_sample_fre);
p=zeros(Num,over_sample_fre);
for n = 1:Num
    d_x = Rot_x-coor(n,1);
    d_y = Rot_y-coor(n,2);
    r = sqrt(d_x.^2+d_y.^2+coor(n,3)^2); % distance
    
    alpha = 1+(1/C_0)*(Rot_x.*d_x+Rot_y.*d_y)./r;% Doppler effector
    
    beta = 1-1i*Rot_v./r./(Fre*2*pi).*((1-alpha-Rot_v^2/C_0^2)...
        ./(Rot_v/C_0*alpha)+1./(Rot_v/C_0*alpha)*(...
        -(Rot_fre*2*pi/C_0^2).*r.*(Rot_x.*d_x+Rot_y.*d_y)...
        /sqrt(Rot_x.^2+Rot_y.^2))); % almost be 0;
    
    d_t(n,:) = floor(r/C_0*over_sample_fre); %time delay;
    p(n,:) = ...
        real(1./(4*pi*r).*exp(1i*2*pi*Fre*t)./alpha.^2.*beta);
        % Signal without delay.
end
%--------------------------------------------------------
% Time delay
P=zeros(size(p));
for i = 1:Num
    for nt = 1:over_sample_fre
        P(i,nt+d_t(i,nt)) = p(i,nt);
    end
end
%--------------------------------------------------------
% Smoothness
for i = 1:Num
    for nt = 1:over_sample_fre
        if P(i,nt) == 0
            P(i,nt) = P(i,nt+1);
        end
    end
end
%******************* Down Sampling ********************************
data = zeros(Num,sample_fre);
for i = 1:sample_fre
    data(:,i) = P(:,(i-1)*over_sample+1);
end
save(data_mat,'data');
clear all;
close all;
clc;
disp('Done data set up');































