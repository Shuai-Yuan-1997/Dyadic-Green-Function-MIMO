% This Matlab code calculates the vector channel matrix of a free-space MIMO communication system based on dyadic Green's function. 
% With this code, it's easy to estimate the capacity or effective degree of freedom (EDOF) of a near-field MIMO communication system.

% Ref: S. S. A. Yuan, Z. He, X. Chen, C. Huang and W. E. I. Sha,"Electromagnetic Effective Degree of Freedom of an MIMO System in Free Space," in IEEE Antennas and Wireless Propagation Letters, vol. 21, no. 3, pp. 446-450, March 2022.

% Input: positions of sources and receivers
% Output: channel matrix considering full polarizations, EDOF

% Shuai S. A. Yuan, Zhejiang University, China
% E-mail: shuaiyuan1997@zju.edu.cn

function vector_channel
clc;
close all;
clear
wavelength=1;  % Notice that the wavelength is set as 1 in this code
wavenumber=2*pi/wavelength;
side_length=10; % Side length of the square source plane and receiving plane
source_num=[2:25];% Number of uniform-distributed sources/receivers along one dimension
for h= 20 % Distance between source and receiving plane
    for k=source_num 

        source_plane=side_length*wavelength; % Side length of source plane
        receiver_plane=side_length*wavelength; % Side length of receiver plane
        distance=h*wavelength; % Distance between source plane and receiver plane
        gap=source_plane/(k); % Gap between ajacent sources
        
        % Number and positions of sources
        N_s=k; % Number of point sources
        x_s=[-(N_s-1)/2:(N_s-1)/2]*gap;
        y_s=[-(N_s-1)/2:(N_s-1)/2]*gap;
        [xx_s,yy_s]=meshgrid(x_s,y_s);
        z_s=-distance/2;
        position_sx=reshape(xx_s,1,N_s^2);
        position_sy=reshape(yy_s,1,N_s^2);
        position_sz=z_s;
        
        % Number and positions of receivers
        N_r=k;% Point receivers number
        x_r=[-(N_r-1)/2:(N_r-1)/2]*gap;
        y_r=[-(N_r-1)/2:(N_r-1)/2]*gap;
        [xx_r,yy_r]=meshgrid(x_r,y_r);
        z_r=distance/2;
        position_rx=reshape(xx_r,1,N_r^2);
        position_ry=reshape(yy_r,1,N_r^2);
        position_rz=z_r;
        
        % Vector channel generation based on the positions of sources and receivers
        for m=1:N_r^2
            for n=1:N_s^2
                channel_matrix(1+(m-1)*3:m*3,1+(n-1)*3:n*3)=dyadic_G(position_rx(m),position_ry(m),position_rz,...
                    position_sx(n),position_sy(n),position_sz);
            end
        end
        correlation=channel_matrix*channel_matrix'; % Correlation matrix
        EDOF=(trace(correlation)/norm(correlation,'fro'))^2; % EDOF calculation
        % Notice that the calculation of EDOF doesn't need normalization. When calculting capacity, proper normalizations should be made! 
        EDOF_matrix(k-1)=EDOF;
    end
    plot(source_num,EDOF_matrix)
    xlabel('number of point sources (along one dimension)')
    ylabel('EDOF')
    hold on
end


function output=dyadic_G(x,y,z,x1,y1,z1) 
% Free-space dyadic Green's function has analytical solution, so it's convenient to model the effect of polarization in near-field MIMO communciations with this method
% x y z denote receiver positions, x1 y1 z1 denote source positions
R=sqrt((x-x1)^2+(y-y1)^2+(z-z1)^2);
wavelength=1;  %
Volume=1;
Z=120*pi;
wavenumber=2*pi/wavelength;
alpha=wavenumber*R;
%  direction
cosx=(x-x1)/R;
cosy=(y-y1)/R;
cosz=(z-z1)/R;
%  costant
const1=1j*wavenumber*Z*wavenumber*Volume*exp(-1j*alpha)/(4*pi*alpha^3);
const2=3-alpha^2+3*1j*alpha;
const3=(alpha)^2-1-1j*alpha;
Greenxx=const1*(const3+cosx*cosx*const2);
Greenyy=const1*(const3+cosy*cosy*const2);
Greenzz=const1*(const3+cosz*cosz*const2);
Greenxy=const1*cosx*cosy*const2;
Greenxz=const1*cosx*cosz*const2;
Greenyz=const1*cosy*cosz*const2;

% Single polarization
% output=[Greenxx 0 0;
%     Greenxy 0 0;
%     Greenxz 0 0]; 

% Full polarizations
output=[Greenxx Greenxy Greenxz;
    Greenxy Greenyy Greenyz;
    Greenxz Greenyz Greenzz];
