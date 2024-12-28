function [ws,wc,wt]=weighting();
% Weight functions WC WS and WT for system
%
%           --------     WS     ------------
%
gainnn=20; % 20
const1= 0.6; 
const2=5.1;
ms=23*gainnn;    
wb1=0.035*const1;  % 0.0001
wb2=0.044*const1;  %  0.0001
wb3=0.27*const1;   %  0.0001
ee=1e-2; 
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
%
%
%           --------     WT     ------------
mt=20*gainnn;    
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


%           --------     WC     ------------
mc=390*gainnn;  %280
wbc=80000; %80000 ;80001 , 90000, ee=1e-2
wbc2=80001;
wbc3=90000;
ee=1e-2; %-4
%ki=0.9;
wc=nd2sys([1 1/mc*wbc],[ee wbc]);
wc2=nd2sys([1 1/mc*wbc2],[ee wbc2]);
wc3=nd2sys([1 1/(mc*1)*wbc3],[ee wbc3]);
%
%wc=daug(wc1,wc2);
wc=daug(wc,wc2,wc3);

%
%**************************************************************************

%**************************************************************************
%                  ----    plot  WS  WC WT    -----
w=logspace(-7,7,100);
%
sv1s=sigma(ws1,w,1);   
sv2s=sigma(ws2,w,1);
sv3s=sigma(ws3,w,1);

sv1t=sigma(wt1,w,1);
sv2t=sigma(wt2,w,1);
sv3t=sigma(wt3,w,1);

sv1c=sigma(wc,w,1);
sv2c=sigma(wc2,w,1);
sv3c=sigma(wc3,w,1);

semilogx(w, 20*log10(sv1s),'b',w, 20*log10(sv2s),'b',w, 20*log10(sv3s),'b' ...
        ,w, 20*log10(sv1t),'m',w, 20*log10(sv2t),'m',w, 20*log10(sv3t),'m', ...
        w, 20*log10(sv1c),'r',w, 20*log10(sv2c),'r',w, 20*log10(sv3c),'r','linewidth',2);
legend('ws','b','c','wt','f','g','wc')
title('WS WT ')
xlabel('rad/sec')
ylabel('dB')
grid
