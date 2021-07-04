function [A,B] = WakeMatrix(N,alpha,beta,Re)
%
% create Orr-Sommerfeld matrices using Chebyshev 
% pseudospectral discretization for plane Poiseuille flow 
% profile 
%
% N  = number of modes
% alpha     = alpha
% beta    = beta 
% Re       = Reynolds number

    global D0 D1 D2 D4

    zi=sqrt(-1);
    
    % mean velocity

    k2 = alpha^2 + beta^2;
    Nos = N+1;
    Nsq = N+1;
    lambda=0.75;
    ymin=-5;ymax=5;
    y=linspace(ymin,ymax,N+1);
    U=1-lambda*exp(-log(2)*y.^2);
    dU=2*lambda*log(2)*y.*exp(-log(2)*y.^2);
    ddU=-2*lambda*log(2)*(2*log(2)*y.^2-1).*exp(-log(2)*y.^2);
        
    % set up Orr-Sommerfeld matrix
    
    B11=D2-k2*D0;
    B12=zeros(Nos);
    A11=alpha*zi*ddU.*D0+(D4-2*k2*D2+k2*k2)/Re-zi*alpha*U.*(D2-k2);
    A12=zeros(Nos);
    er  = -200*zi;
    A11 = [er*D0(1,:); er*D1(1,:); A11(3:Nos-2,:); er*D1(Nos,:); er*D0(Nos,:) ];
    B11 = [D0(1,:); D1(1,:); B11(3:Nos-2,:); D1(Nos,:); D0(Nos,:)];  
       
    % set up Squire matrix and (cross-term) coupling matrix
    
    A21=-zi*beta*dU.*D0;
    A22=(D2-k2)/Re-zi*alpha*U.*ones(1,length(U)).*D0(1:Nos,:);
    B21=zeros(Nsq);
    B22=D0;
    A22 = [er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
    A21 = [zeros(1,Nos); A21(2:Nsq-1,:); zeros(1,Nos)];
    
    % combine all the blocks 
    A = [A11 A12; A21 A22];
    B = [B11 B12; B21 B22];
    