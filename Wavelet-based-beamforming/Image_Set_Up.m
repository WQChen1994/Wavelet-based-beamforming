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
%********************** Imaging Set Up ********************************
clear all;close all;clc;
C_0 = 340;              % Speed of sound.
sample_fre = 40000;     % Sampling frequency Hz
Fre = 2000;             % signal frequency
Num = 62;               % Number of sensors 
Rot_Fre = 50;
%----------------------------------------------------
% load data
main_root = 'D:\MATLAB\CWT_fust\Wavelet_beamforming';
CWT_path = '\CWT';
CWT_mat = [main_root,CWT_path,'\CWT_rotating.mat'];
Fre_mat = [main_root,'\F_cwt.mat'];
% sensor coordinates
mics_coor_mat = ...
    [main_root,'\mics_coor_Array.mat'];
load(CWT_mat); 
load(mics_coor_mat);          % Array coordinates --> coor
load(Fre_mat);
disp('Done load');
%----------------------------------------------
T_begin = floor(0.3*sample_fre);          % Begining Imaging time.
T_during = floor(0.1*sample_fre);         % Beamforming time.
T_end=T_begin+T_during;                   % Stop Imaging time.
%-----------------------------------------------
% Beamforming Set Up
r_begin = 0;                % begin imaging location
r_end = 0.5;                   % end imaging location
theta_begin = 0;
theta_end =2*pi;
Num_r = 50;                    % Number of x mesh
Num_t = 50;
ra = linspace(r_begin,r_end,Num_r);
theta = linspace(theta_begin,theta_end,Num_t);
%-----------------------------------------------
for nr = 1:Num_r
    for nth = 1:Num_t
        
        d_x(:) = ra(nr)*cos(theta(nth))-coor(:,1);
        d_y(:) = ra(nr)*sin(theta(nth))-coor(:,2);  % distance
        dr(:,nr,nth) = sqrt(d_x.^2+d_y.^2+coor(:,3)'.^2);% distance
        d_ux = ra(nr)*sin(theta(nth))*2*pi*Rot_Fre;
        d_uy = -ra(nr)*cos(theta(nth))*2*pi*Rot_Fre;
        v = ra(nr)*2*pi*Rot_Fre;    % velocity
        dt(:,nr,nth) = floor(dr(:,nr,nth)./C_0*sample_fre)+1;
        alpha(:,nr,nth) = 1+(d_ux.*d_x+d_uy.*d_y)./dr(:,nr,nth)'/C_0;
        fre_m = Fre./alpha(:,nr,nth);
        for n = 1:Num
            [~,fre_x(n,nr,nth)] = min(abs(F_cwt-fre_m(n)));
        end
        clc;
        disp(['Seting   ',num2str(nr*Num_t-Num_t+nth),...
            '   /   ',num2str(Num_r*Num_t)]);
    end
end

%**********************************************************************
%************************** Imaging    ********************************
disp('Imaging');
SPL = zeros(Num_r,Num_t,T_during);
for nt =1:T_during
    tic
    for nr = 1:Num_r
        for nth = 1:Num_t
            for n = 1:Num
                Ya(n,1) = P_cwt(fre_x(n,nr,nth),nt+dt(n,nr,nth),n);
                % Delay time and Choose frequency
            end
            % Following Classical beamforming
            G = 1./(4*pi*dr(:,nr,nth))./alpha(:,nr,nth).^2;
            CSM = Ya*Ya';
            w = G'/norm(G);
            SPL(nr,nth,nt) = 10*log10(abs(w*CSM*w')/4*10^10);
        end
    end
    clc
    disp('Imaging');
    disp([' Done beamforming from time  ',num2str(T_begin/sample_fre),...
        '  to time ', num2str(T_end/sample_fre),...
        '  ',num2str(nt),'  /   ',num2str(T_during)]);
    toc
end
Result_path='\Result';
Result_mat=[main_root,Result_path,'\result_rotating.mat'];
save(Result_mat,'SPL','-v7.3')
%--------------------------------------------------------
% Show
[Ra,Theta]=meshgrid(ra,theta);
xa=Ra.*cos(Theta);
ya=Ra.*sin(Theta);              % mesh polar-->Decare
for nt=1:T_during
    surf(xa,ya,SPL(:,:,nt)');shading interp;view(2)
    caxis([max(max(SPL(:,:,nt)))-12,max(max(SPL(:,:,nt))) ])
    colorbar;
    colormap jet;
    title([num2str(Fre),' Hz   ',...
        'Time=   ',num2str(nt),'   /   ',num2str(T_during)]);
    drawnow;
end




















