%**************************************************************
%System identification with least square and Kalman Filting methods
%...........................................................
% Written  by Wissem CHIHA 
%Tunisia Polytechnic School 
%wissem.chiha@ept.ucar.tn
%***********************************************************

%-------------start script---------------------------------

%system a identifier: first order sytem : H(s)=k/(1+tau*s)
%method : least square method 
%simulate system
% Define the transfer function
clc
clear 
k = 1;           % Static Gain
tau = 1.4;         % Time constant
syms s           % laplace variable 
H = k/(1+tau*s); % Transfer function
% Plot the step response
uo=1;
T=30;            % simulation time duration 
syms t           % time variable 
v=ilaplace(H*(uo/s),t);

% Define the length of the white noise signal
n = 800;
% Generate white noise using randn
white_noise = randn(n,1);
noise_variance=var(white_noise);
noise_mean_value=mean(white_noise);
noise_mangitude=0.04;
% generate the noised output signal
Ts=T/(n-1);   %sampling  period  
t=0:Ts:T;
y=subs(v)'- noise_mangitude.*white_noise;

figure(1)
subplot(211)
plot(t,subs(v)','r',t,y),grid 
title({'First order system with parameters',
        ['k = ',num2str( k),', tau= ',num2str(tau)]})
legend('True output','Noised output');   
subplot(212)
plot(t,y-subs(v)','g'),grid
title( {' Mesurement Error   ';
        ['noise mean = ',num2str( noise_mean_value),', noise variance= ',num2str( noise_variance)]})
    
%**************************************************************
%...................Recursive least squares method.............
%**************************************************************
u=uo.*ones(n,1);
theta_init=[1,0.3];
% Parameter initialization
theta = theta_init;
P = diag([1, 1]);
y_est=0;
N=100;        %  mesures sample time is N*Ts 
y=y(1:N:n);
% Loop to update parameters with each new observation
for k = 2:length(y)

    % Calculation of the system estimate at time k
    y_est = exp(-Ts/theta(2)) * y_est + (1 - exp(-Ts/theta(2))) * theta(1) * u(k-1);
    %Calculation of prediction error
    e = y(k) - y_est;
    % Calculation of update matrices
    K = P * [1; u(k-1)] / (1 + [1; u(k-1)]' * P * [1; u(k-1)]);
    P = (eye(2) - K * [1, u(k-1)]) * P;
    % Parameter update
    theta = theta + K * e;
end
theta=double(theta);
%plot the response of the new  estimated sytem and error
figure(2)
subplot(211)
syms p
sym=theta(1,2)/(theta(1,1)*p +1);
y_n=ilaplace(sym*uo/p);
plot(t,subs(y_n)','r',t,subs(v)','b'),grid
title({'parameters estimation with least recursive square method ';
        ['k = ',num2str( theta(1,2)),', tau= ',num2str(theta(1,1))]}) 
legend('True output','Estimeted output');
subplot(212)
plot(t,subs(y_n)'-subs(v)','c'),grid
title('Estimation Error RLS')

%Display results 
disp('The estimated parameterswith least recursive square method are:');

disp(double(theta));

%******************************************************************
%.........................Kalman filtring method..................
%******************************************************************

%Define the state space form of the system 
A=-tau;
B=k;
C=1 ;
D=0 ;
Ts=-1;
%create the sytem dynamics with the input noise w 
sys = ss(A,[B B],C,D,Ts,'InputName',{'u' 'w'},'OutputName','yt'); 
%Define the process noise covariance and the sensor noise covariance 
Q = 0.5; 
R = 0.5; 
%we use kalman Matlab command to design the filter 
[kalmf,L,~,Mx,Z] = kalman(sys,Q,R);
%we eliminate the state estimation form thee filter 
kalmf=kalmf(1,:);

%..............testing the filter on th system data.................
vIn = sumblk('y=yt+v');
kalmf.InputName = {'u','y'};
kalmf.OutputName = 'ye';
SimModel = connect(sys,vIn,kalmf,{'u','w','v'},{'yt','ye'});
%generate the step input and noises vectors 
t=(1:T)';
u=uo.*ones(length(t),1);
w = sqrt(Q)*randn(length(t),1);
v = sqrt(R)*randn(length(t),1);
%simulate the model and get the results 
out = lsim(SimModel,[u,w,v]);
%extraction and plotting results 
yt=out(:,1) ;      %the true response of our  system
ye=out(:,2) ;      %the kalman filtred respone
y=yt+v;           %the noised response 

figure(3)
subplot(211), plot(t,yt,'b',t,ye,'r--'), 
xlabel('Number of Samples'), ylabel('Output')
title('Kalman Filter Response')
legend('True output','Filtered output')
subplot(212), plot(t,yt-y,'g',t,yt-ye,'r--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered')



