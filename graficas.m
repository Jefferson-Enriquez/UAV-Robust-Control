close all
w=logspace(-4,4,10000);
sensi=pck(al-bl*cl, bl, -cl, eye(3));
tsensi=pck(al-bl*cl, bl, cl, 0*eye(3));
csensi=pck(al-bl*cl, bl, cl, 0*eye(3));
%
sv11 = sigma(ss(al-bl*cl, bl, -cl, eye(3)),w);
sv11 = 20*log10(sv11);
figure;
semilogx(w, sv11,'--','linewidth',2)
%% Sensitivity S
title('x,y and z speeds')
grid
xlabel('Frequency (radian/second)')
ylabel('Amplitude (decibel)')
hold on


%           --------     WS     ------------
%
gainnn=20; % 20
const1= 0.6; % 5.5 para drone con valores del paper
const2=5.1;
ms=23*gainnn;    % 45  3.3 garante um overshot Mp < 6dB
wb1=0.035*const1;  % 0.0001
wb2=0.044*const1;  %  0.0001
wb3=0.27*const1;   %  0.0001
ee=1e-2; %10   6  con 2 da un mejor valor pero es bueno?
ki=1;
num1=conv([1/sqrt(ms) 10^(ki-1)*wb1],[1/sqrt(ms) 10^(ki-1)*wb1]);
den1=conv([1 10^(ki-1)*wb1*sqrt(ee)],[1 10^(ki-1)*wb1*sqrt(ee)]);
ws1=nd2sys(num1,den1);

num2=conv([1/sqrt(ms) 10^(ki-1)*wb2],[1/sqrt(ms) 10^(ki-1)*wb2]);
den2=conv([1 10^(ki-1)*wb2*sqrt(ee)],[1 10^(ki-1)*wb2*sqrt(ee)]);
ws2=nd2sys(num2,den2);

num3=conv([1/sqrt(ms) 10^(ki-1)*wb3],[1/sqrt(ms) 10^(ki-1)*wb3]);
den3=conv([1 10^(ki-1)*wb3*sqrt(ee)],[1 10^(ki-1)*wb3*sqrt(ee)]);
ws3=nd2sys(num3,den3);
%
ws=daug(ws1,ws2,ws3);

sv1s=sigma(ws1,w,1);
sv2s=sigma(ws2,w,1);
sv3s=sigma(ws3,w,1);
semilogx( w, (20*log10(sv1s)),w, 20*log10(sv2s),w, 20*log10(sv3s),'linewidth',2);
legend('Sx','Sy','Sz','WSx','WSy','WSz')





% -----------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------


sensi=pck(al-bl*cl, bl, -cl, eye(3));
tsensi=pck(al-bl*cl, bl, cl, 0*eye(3));
csensi=pck(al-bl*cl, bl, cl, 0*eye(3));

sv = sigma(ss(al-bl*cl, bl, cl, 0*eye(3)),w);
sv = 20*log10(sv);
figure;
semilogx(w, sv,'--','linewidth',2)
%% Complementary Sensitivity T
title('x,y and z speeds' )
grid
xlabel('Frequency (radian/second)')
ylabel('Amplitude (decibel)')
hold on
%
%           --------     WT     ------------
mt=20*gainnn;    % 15 1.28garante um overshot Mp < 2dB
wbt=900*const2; % 400*const
wbt2=900*const2; % 400*const
wbt3=1000*const2+6; % 2.4
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
semilogx( w, (20*log10(sv1t)),'b',w, 20*log10(sv2t),w, 20*log10(sv3t),'linewidth',2);
legend('Tx','Ty','Tz','WTx','WTy','WTz')


% -----------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%------- Plot the singular Values WsS 
wssensi=mmult(ws,sensi);
wtsensi=mmult(wt,tsensi);
wcsensi=mmult(wc,csensi);

sv1= sigma(wssensi,w);
figure;
semilogx(w,sv1,'linewidth',2)

%title('Ws*S' )
title('||W_s S|| < 1' )
grid
xlabel('Frequency (radian/second)')
ylabel('Magnitude')
legend('V_x','V_y','V_z')