clear;
close all;
clc;

xx=[0.01,0.01,0.0001,0.01,0.01,0.01,0.01,0.0001,0.01,0.01];
uu=[0.00001,0.00001,0.00001];
%   Vx, Vy, Vz, 
int5=[0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001] ;%int5=[0,0,0,0.1,0.1,0.1,0.1,0,0.1,0.1] ;

[A,B,C,D]=linmod('PlantUAV',xx,uu);
   [A1,B1,C1,D1]=linmod('PlantUAV',xx,uu);

AAA=A;  
A=A(4:10,4:10);
BBB=B;
B=B(4:10,:);
CCC=C;
C=[0 0 0 0 0 1 0 0 0 0;
   0 0 0 0 0 0 1 0 0 0;
   0 0 0 0 0 0 0 1 0 0]
C=C(:,4:10)
% %---% reduction of B
% %---% reduction of D 
D=zeros(3,3);

%observability and controllability
OBB=obsv(A,C);
COO=ctrb(A,B);
obsv_original=rank(OBB)
crt_original=rank(COO)


% Pack G and Weightings

[ws,wc,wt]=weighting();


hold on
ww=logspace(-7,7,100);
sv = sigma(ss(A,B,C,D),ww);
sv = 20*log10(sv);
semilogx(ww, sv,'linewidth',2)
xlabel('rad/sec')
ylabel('dB')


% Augment Plant 

[ns nc] = size(B);  % ns = number of inputs;  nc = number of controls;   
var=0;
aa = [ A            B
      var*ones(nc,ns)    var*ones(nc,nc) ];
bb = [ var*ones(ns,nc)
      eye(nc)      ];
cc = [ C  var*ones(nc,nc) ];
dd = var*ones(nc,nc);

OBB=obsv(aa,cc);
COO=ctrb(aa,bb);
obs_ampli=rank(OBB)
crt_ampli=rank(COO)


%G=pck(aa,bb,cc,dd);  %expanded plant
G=pck(A,B,C,D);       %original plant
systemnames = 'G wc ws wt';
%3 inputs and 3 outputs are being considered
inputvar='[r(3);u(3)]'; % linear speed control
outputvar='[ws;wc;wt;r-G]';

input_to_G='[u]';
input_to_wc='[u]';  
input_to_ws='[r-G]';
input_to_wt='[G]';
sysoutname='pext';            % expanded plant
cleanupsysic='yes';
sysic;


% synthesis of the controller

%
	nmeas=3;        % number of measured outputs
	ncon=3;         % number of control inputs
	gmin=0.0001;
	gmax=10000000;
	tol=0.001;

[a_pext,b_pext,c_pext,d_pext]=unpck(pext);
sys_pext=ss(a_pext,b_pext,c_pext,d_pext);
P1=minreal(sys_pext);
[K1,Nsc1,gamma(1),info]=hinfsyn(P1,nmeas,ncon,'method','ric','Tolgam',1e-3);

% [K1,Nsc1,gamma(1),info]=hinfsyn(P1,nmeas,ncon,'method','ric','Tolgam',1e-3)
 
sysk=K1;
[ak1,bk1,ck1,dk1]=ssdata(sysk);

KA=ak1;
KB=bk1;
KC=ck1;
KD=dk1;

%-----------------Result analysis
k=pck(ak1,bk1,ck1,dk1);
tzw=Nsc1;
gsuplot=gamma;
fprintf('The gamma suboptimal = %4.5f\n',gsuplot)
sysh=mmult(G,k);
[al,bl,cl,dl]=unpck(sysh);
[ak1,bk1,ck1,dk1]=unpck(k);

%-------- Closed loop analysis 
clpoles=eig(al-bl*cl) ; % Poles in closed loop
clsys=ss(al-bl*cl, bl, cl, 0*eye(3)); % closed loop system

figure;
pzmap(clsys);
grid on
clzeros=zero(clsys);
clpoles=pole(clsys);

% To see if there are zeros close to or superpositioned over poles
clzerosStr=num2str(clzeros);
clpolesStr=num2str(clpoles);
disp('see the pzmap and indentify the location of poles and zeros')
disp('and press Return to continue')
%pause
%
w=logspace(-7,7,100);
sensi=pck(al-bl*cl, bl, -cl, eye(3));
tsensi=pck(al-bl*cl, bl, cl, 0*eye(3));
csensi=pck(al-bl*cl, bl, cl, 0*eye(3));
%
sv11 = sigma(ss(al-bl*cl, bl, -cl, eye(3)),w);
sv11 = 20*log10(sv11);
figure;
semilogx(w, sv11)
%clear sv
title('Sensitivity')
grid
xlabel('Frequency (rad/s)')
ylabel('Singular Values (dB)')
hold on

%           --------     WT     ------------
const=2;
mt=20;    % 15 1.28garante um overshot Mp < 2dB
wbt=900*const; % 400*const
wbt2=900*const; % 400*const
wbt3=1000*const+6; % 2.4
ee=1e-4; %-4
num1=conv([1 wbt/sqrt(mt)],[1 wbt/sqrt(mt)]);
den1=conv([sqrt(ee) wbt],[sqrt(ee) wbt]);
%wt1=nd2sys([1 1/mt*wbt],[ee wbt]);
num2=conv([1 wbt2/sqrt(mt)],[1 wbt2/sqrt(mt)]);
den2=conv([sqrt(ee) wbt2],[sqrt(ee) wbt2]);
num3=conv([1 wbt3/sqrt(mt)],[1 wbt3/sqrt(mt)]);
den3=conv([sqrt(ee) wbt3],[sqrt(ee) wbt3]);

wt1=nd2sys(num1,den1);
wt2=nd2sys(num2,den2);
wt3=nd2sys(num3,den3);
wt=daug(wt1,wt2,wt3);

sv1t=sigma(wt1,w,1);
sv2t=sigma(wt2,w,1);
sv3t=sigma(wt3,w,1);
semilogx( w, (20*log10(sv1t)),'m',w, 20*log10(sv2t),'m',w, 20*log10(sv3t),'m','linewidth',2);


