% function[xdot]=lineal44simu(x,Vx,y,Vy,z,Vz,phi_l,phi_d_l,tetha_l,tetha_d_l,Fx,Fy,Fz)
function[xdot]=lineal5simu(x,y,z,phi_l,tetha_l,Vx,Vy,Vz,phi_d_l,tetha_d_l,Fx,Fy,Fz)
%-------------------------MATRICES---------------- 
%constantes
ml=0.2; mc=2.5; L=2; g=9.8;  % datos reales del paper
% ml=0.0001; mc=0.7; L=0.0001; g=9.8; %considerando solo dron , estabiliza en 2 s
%ml=0.1; mc=0.7; L=0.1; g=9.8; %considerando dron con pendulo libiano da
% da un gama de 13
yaw=0; % no es una entrada es cosntante

nu=[Vx; Vy; Vz; phi_d_l; tetha_d_l] ;

%%%%%%    M_n  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=0.001;
M_n=[ml+mc 0 0 0 L*ml*cos(tetha_l);
     0 ml+mc 0 -L*ml*cos(phi_l)*cos(tetha_l) L*ml*sin(phi_l)*sin(tetha_l);
     0 0 ml+mc -L*ml*cos(tetha_l)*sin(phi_l) -L*ml*cos(phi_l)*sin(tetha_l);
     0 -L*ml*cos(phi_l)*cos(tetha_l) -L*ml*cos(tetha_l)*sin(phi_l) (L^2)*ml*(cos(tetha_l))^2+E*(sin(tetha_l))^2 0;
     L*ml*cos(tetha_l) L*ml*sin(phi_l)*sin(tetha_l) -L*ml*cos(phi_l)*sin(tetha_l) 0 (L^2)*ml];

%%%%%%    C_nv  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_nv=[0 0 0 0 -L*tetha_d_l*ml*sin(tetha_l);
      0 0 0 L*ml*(phi_d_l*cos(tetha_l)*sin(phi_l)+tetha_d_l*cos(phi_l)*sin(tetha_l)) L*ml*(phi_d_l*cos(phi_l)*sin(tetha_l)+tetha_d_l*cos(tetha_l)*sin(phi_l));
      0 0 0 -L*ml*(phi_d_l*cos(phi_l)*cos(tetha_l)-tetha_d_l*sin(phi_l)*sin(tetha_l)) -L*ml*(tetha_d_l*cos(phi_l)*cos(tetha_l)-phi_d_l*sin(phi_l)*sin(tetha_l));
      0 0 0 -0.5*(L^2)*tetha_d_l*ml*sin(2*tetha_l)+0.5*E*tetha_d_l*sin(2*tetha_l) -0.5*(L^2)*phi_d_l*ml*sin(2*tetha_l);
      0 0 0 0.5*(L^2)*phi_d_l*ml*sin(2*tetha_l) 0];

%%%%%%    G_n  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_n=[0;
    0;
    -g*(ml+mc);
    L*g*ml*cos(tetha_l)*sin(phi_l);
    L*g*ml*cos(phi_l)*sin(tetha_l)];

%%%
%%%    D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=2; % 1.2
D=diag([0 0 0 d d]);
% %%%%%%   R_dot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  R=[cos(yaw)*cos(pitch)-sin(roll)*sin(yaw)*sin(pitch), -cos(roll)*sin(yaw), cos(yaw)*sin(pitch)+cos(pitch)*sin(roll)*sin(yaw) ;...
%         cos(pitch)*sin(yaw)+cos(yaw)*sin(roll)*sin(pitch),  cos(roll)*cos(yaw), sin(yaw)*sin(pitch)-cos(yaw)*cos(pitch)*sin(roll);...
%         -cos(roll)*sin(pitch), sin(roll), cos(pitch)*cos(roll)]   ;  % ecuacion 33
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%    tau  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=[0,0,fz_in]';
% F=R*f;
% tau=[F' 0 0]';
tau=[Fx Fy Fz 0 0]';

%%%%%%    tau_a  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tau_a=[wr;wr;wr;0;0];
tau_a=[0;0;0;0;0];

%PARA LAS VARIABLES DE ESTADO:
Nu_dot= inv(M_n)*(-C_nv*nu-G_n-D*nu+tau+tau_a); % 
%-------------------------------------------------
% x1=x; x2=Vx; x3=y; x4=Vy; x5=z; x6=Vz;
% x7=phi_l; x8=phi_d_l; x9=tetha_l; x10=tetha_d_l;
%varibles de estado finales
% xdot= [Vx;
%        Nu_dot(1);
%        Vy;
%        Nu_dot(2);
%        Vz;
%        Nu_dot(3);
%        phi_d_l;
%        Nu_dot(4);
%        tetha_d_l;
%        Nu_dot(5)
%        ];
xdot=[nu;
      Nu_dot];