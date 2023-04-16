 
%*****************************************************************************
%%  PARAMETRIC SYSTEM IDENTIFICATION 
%
%  Discussion:
%
%    This program assumes that you have a Controlled Auto Regressive
%    model of your system described by the form :
% 
%        A(z)*Y=B*U+E
%   
%  Usage:
%
%        identification ( 'data', na, nb )
%
%  where:
%        data.mat :Experimental response values of the system 
%        na       : denominator Order 
%        nb       :Numerator order 
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Author:
%
%     Wisssem Chiha 
%     wissem.chiha@ept.ucar.tn
%
%%
clc
clear all 
 %%
 % SET PROGRAM GENERL PRAMETERS:
       x=[-1.67 1.33 -0.07 -1.21 2.44 1.92 0.38 -0.2 1.0 0.89 ];
       y=[-0.48 1.15 0.31 -1.01 2.2 2.61 1.75 0.86 0.76 1.86];
       na=2;  nb=1;
       lamda=0.7;
 
       iter_number=5;  %iteration number for GLS
       q=nb;           %filter order 

%%
%**************************************************************
%................Simple Least squares method...................
%**************************************************************

if na >= nb+1

N=size(x,2); 
if N > min(na,nb+1)     

H=zeros(N,na+nb+1);
   for i=1:N-na+1
       for j=1:na
           H(i,j)=-y(na-j+i);
       end
   end
 for i=1:N-na+1
     for j=1:nb+1
          H(i,j+na)=x(na-j+i);
     end
 end
theta=inv(H'*H)*((H')*y');
disp('**********************************************');
disp('Simple Least squares Estimation Parameters: ');
disp('**********************************************');
disp(['A(q)=[a0...a',num2str(na),']']);
disp(theta(1:na)');
disp(['B(q)=[b0...b',num2str(nb),']']);
disp(theta(na+1:na+nb+1)');
else
    error('Error lack of Samples values ');
end
else
    error('Invalid System Model or Model not Supported : non causal system detected :try na > nb');
     
     
end 
%%
%**************************************************************
%................Weighted Least squares method...................
%**************************************************************
lamda=0.7;
Q=zeros(N);
for i=1:N
    Q(i,i)=lamda^(i-1);
end 
theta=inv(H'*H)*((H')*(Q*y'));
disp('**********************************************');
disp('Weighted Least squares  Estimation Parameters:');
disp('**********************************************');
disp(['Weight factor: ',num2str(lamda)]);

disp(['A(q)=[a0...a',num2str(na),']']);
disp(theta(1:na)');
disp(['B(q)=[b0...b',num2str(nb),']']);
disp(theta(na+1:na+nb+1)');

%%
%%**************************************************************
%................Recursive Least squares method.................
%***************************************************************
%forgetting factors
lamda=0.99;
%Regularization parameter
alpha=1e-6;
%Parameters Initilisation
P=alpha.*eye(na+nb+1);
theta_RLS1=zeros(na+nb+1,1);
%Main Loop
for k=1:N-na+1
    yk=y(k);
    
   % Construction of the regression vector h
    hk= zeros(na+nb+1, 1);
    for i = 1:na
        hk(i, 1) = -y(na-i+k);
    end
    for i = 1:nb+1
        hk(na+i, 1) = x(na-i+k);
    end
    %Update paramters
    G = (P * hk) / (lamda +  hk' * P * hk);
    P = (P - G * hk' * P)./lamda ;
    theta_RLS1 = theta_RLS1 + G * (yk - hk' * theta_RLS1);
end
disp('**********************************************');
disp('Recursive Least squares  Estimation Parameters: ');
disp('**********************************************');
disp(' ');
disp(['Algorithm with forgetting factor :',num2str(lamda),]);

disp(['A(q)=[a0...a',num2str(na),']']);
disp(theta_RLS1(1:na)');
disp(['B(q)=[b0...b',num2str(nb),']']);
disp(theta_RLS1(na+1:na+nb+1)');

%check  the P matrix not explosed 
%Calculate the spectral radius 
spectral_radius = max(abs(eig(P)));
if spectral_radius > 1
   warning(['The spectral radius of matrix P is ', num2str(spectral_radius)]); 
   warning('Itervative process might be  divergent !');
end
%%
%-------------------------------------------------------------------------
%forgetting factors
alpha_tr_cte=0.1; %(=lamda1/lamda2)
tau=1; %gain 
%Initialisation
P=alpha.*eye(na+nb+1);
theta_RLS2=zeros(na+nb+1,1);
for k=1:N-na+1
     yk=y(k);
  
    % Construction of the regression vector h
    hk= zeros(na+nb+1, 1);
    for i = 1:na
        hk(i, 1) = -y(na-i+k);
    end
    for i = 1:nb+1
        hk(na+i, 1) = x(na-i+k);
    end
    lamda1=(1/(tau*N))*trace(P-( P *(hk * hk') * P)*inv(alpha_tr_cte +hk' * P *hk));
    G=(P * hk)*inv(alpha_tr_cte+hk' * P * hk);
    theta_RLS2=theta_RLS2+ G *(yk-hk' * theta_RLS2);
    P=1/lamda1*(eye(size(G*hk'))-G * hk') * P;
end
 disp(['Algorithm with constant  trace : alpha =',num2str(alpha_tr_cte),]);
 disp(['A(q)=[a0...a',num2str(na),']']);
disp(theta_RLS2(1:na)');
disp(['B(q)=[b0...b',num2str(nb),']']);
disp(theta_RLS2(na+1:na+nb+1)');

%%
%%**************************************************************
%................General Least squares method.................
%***************************************************************
 

 for i=1:iter_number
   %First MCO estimation 
   H=zeros(N,na+nb+1);
   for i=1:N-na+1
       for j=1:na
           H(i,j)=-y(na-j+i);
       end
   end
 for i=1:N-na+1
     for j=1:nb+1
          H(i,j+na)=x(na-j+i);
     end
 end
theta_k=inv(H'*H)*((H')*y');
%--------------------------------------------------
%compute Eg_k
Eg_k=y'-H*theta_k;
disp(size((Eg_k)));
%--------------------------------------------------
%Filetr paramter estimation 
Hf=zeros(nb,q);
for t=1:nb
    for j=1:q
        Hf(t,j)=-Eg_k(nb+na-t-j+1);
    end
end
theta_f=inv(Hf'*Hf)*(Hf'*Eg_k);
%-------------------------------------------------
%filtrage
for k=q+1:N
 for s=1:q
     y(k)=y(k)+theta_f(s)*y(k-s);
     x(k)=x(k)+theta_f(s)*x(k-s);
end
end
 end
disp('***************************************************');
disp('General Least square method estimation parameters :');
disp('***************************************************');
disp(' ');
disp(['fileter order :',num2str(q)])
disp(['iteration number : ',num2str(iter_number)]);
disp(' ');
disp(['A(q)=[a0...a',num2str(na),']']);
disp(theta_k(1:na));
disp(['B(q)=[b0...b',num2str(nb),']']);
disp(theta_k(na+1:na+nb+1)');




 
 
 
