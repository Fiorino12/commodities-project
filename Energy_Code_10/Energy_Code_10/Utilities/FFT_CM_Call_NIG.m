function [Price]=FFT_CM_Call_NIG(Strike,T,r,S0,CharFunc)
% Price of a Plain Vanilla Call 
% exploiting the Carr-Madan algorithm

% discretization parameter
%Npow = 14; N=2^Npow; A = 600; 
%Npow=16; N=2^Npow; A=1000;
Npow = 16; N=2^Npow; A = 1000;

% v-> compute integral as a summation
eta=A/N; v=[0:eta:A*(N-1)/N]; v(1)=1e-22;
% lambda-> compute summation via FFT
lambda=2*pi/(N*eta); 
k=-lambda*N/2+lambda*(0:N-1);

% Fourier transform of z_k
% CharFunc=@(v) exp(T*CharExp(v,params));
% disp('RiskNeutral Check')
% CharFunc(-1i)
Z_k=(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));  %exp(1i*r*v*T).*...
% Option Price
w=ones(1,N); w(1)=0.5; w(end)=0.5;
x=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(x)/pi);
%C=S0*(z_k+max(1-exp(k-r*T),0));
C=exp(-r*T)*S0*(z_k+max(1-exp(k),0));
K=S0*exp(k);
% Output
index=find( K>0.1*S0 & K<3*S0 );
%length(index)
C=C(index); K=K(index);
% plot(K,C)
% title( 'Option Price' );
% xlabel('Strike');
Price=interp1(K,C,Strike,'spline');
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% function V=CharExp(v,params)
% % risk-neutral characteristic exponent
% sigma=params(1);
% theta=params(2);
% k=params(3);
% V=@(u) -log(1+u.^2*sigma^2*k/2-1i*theta*k*u)/k;
% drift_rn=-V(-1i); % Drift Risk_neutral
% V=drift_rn*1i*v+V(v);















