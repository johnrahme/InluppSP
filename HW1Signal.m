% Signalbehandling HW 1 %
% Vocal - near end signal
% Drum  - far  end signal
clear all
close all

[V,Fs_V] = audioread('vocal.wav');
[DL,Fs_D] = audioread('drumloop.wav');

phi = zeros(length(V),1);
theta = zeros(length(V),1);
delay = [0.05 0.1 0.3].*Fs_D;
c=[0.1 0.4 0.2];

%Same length
DL(length(DL)+1:length(V))=0;

%construct echo
for n = delay(1) + 1:length(DL)+delay(1)
    y1(n)=c(1)*DL(n-delay(1));
end

for n = delay(2) + 1:length(DL)+delay(2)
    y2(n)=c(2)*DL(n-delay(2));
end

for n = delay(3) + 1:length(DL)+delay(3)
    y3(n)=c(3)*DL(n-delay(3));
end
%sum of echos
y1(length(y1)+1:length(y3)) = 0;
y2(length(y2)+1:length(y3)) = 0;
y=(y1+y2+y3)';

t1=1:length(y);
t2=1:length(DL);
figure;
subplot(2,1,1);
plot(t1,y);
subplot(2,1,2);
plot(t2,DL);

N = sqrt(0.003)*randn(length(V),1);
U = V+N;
t_u = 1:length(U);
figure(2);
subplot(2,1,1);
plot(t_u,U)
subplot(2,1,2);
plot(t_u,V);

%Sum of noise and echo
U(length(U)+1:length(y)) = 0;
s = y+U;
t_s = 1:length(s);
figure(3);
plot(t_s,s);



%-------------------------------------------------------------------------
%Estimation of c usin LMS

c_hat = zeros(3,1);
my = 0.01;
DL(length(DL)+1:length(s)) = 0;
y_hat=zeros(length(DL),1);
e=zeros(length(DL),1);
phi=zeros(3,1);
for n = 1:length(DL)
    for m = 1:3
        if n-delay(m)>0
            phi(m)=DL(n-delay(m));
        else
            phi(m)=0;
        end
    end
    y_hat(n)=transpose(c_hat)*phi;
    e(n)=s(n)-y_hat(n);
    c_hat=c_hat+2*my*phi*e(n);
end

disp('LMS known: ')
disp(c_hat);
%Estimation of c using RLS

c_hatRLS1 = zeros(3,1);
%A lot of disturbance since white noise, we need a small rho
rho = 0.01;
N = 3;
P = rho*eye(N);
lambda = 0.98;
phi1 = zeros(3,1);
e1 = zeros(length(DL),1);
for n = 1:length(s)
    for m = 1:3
        if n-delay(m)>0
            phi1(m) = DL(n-delay(m));
        else
            phi1(m) = 0;
        end
    end
    
    %P(n) = (1/lambda)*(P(n)-(P(n)*phi*(phi')*P(n))/(lambda+((phi')*P(n)*phi)));
    %P=(P-(P*phi*transpose(phi)*P)/(lambda+transpose(phi)*P*phi))/lambda
    %K = P*phi;
    K1 = P*phi1/(lambda+(phi1')*P*phi1);
    e1(n) = s(n)-c_hatRLS1'*phi1;
    c_hatRLS1 = c_hatRLS1+transpose(K1'*e1(n));
    
end
disp('RLS known: ')
disp(c_hatRLS1)


%-------------------------------------------------------------------------

% LMS on unknown delay..
N_M = 5;
delay1 = (1/Fs_D:5000/Fs_D:N_M)*Fs_D;
c_hatLMS = zeros(length(delay1),1);
my = 0.01;
DL(length(DL)+1:length(s)) = 0;
y_hat=zeros(length(DL),1);
e=zeros(length(DL),1);
phi=zeros(length(delay1),1);
for n = 1:length(DL)
    for m = 1:length(delay1)
        if floor(n-delay1(m))>0
            phi(m)=DL(floor(n-delay1(m)));
        else
            phi(m)=0;
        end
    end
    y_hat(n)=transpose(c_hatLMS)*phi;
    e(n)=s(n)-y_hat(n);
    c_hatLMS=c_hatLMS+2*my*phi*e(n);
end

disp('LMS unknown');
disp(c_hatLMS)

% RLS on unknown delay

c_hatRLS = zeros(length(delay1),1);
N = length(delay1);
P = rho*eye(N);
phi = zeros(length(delay1),1);
e = zeros(length(DL),1);
for n = 1:length(s)
    for m = 1:length(delay1)
        if floor(n-delay1(m))>0
            phi(m)=DL(floor(n-delay1(m)));
        else
            phi(m) = 0;
        end
    end
    
    %P(n) = (1/lambda)*(P(n)-(P(n)*phi*(phi')*P(n))/(lambda+((phi')*P(n)*phi)));
    %P=(P-(P*phi*transpose(phi)*P)/(lambda+transpose(phi)*P*phi))/lambda
    %K = P*phi;
    K = P*phi/(lambda+(phi')*P*phi);
    e(n) = s(n)-c_hatRLS'*phi;
    c_hatRLS = c_hatRLS+transpose(K'*e(n));
    
end

disp('RLS unknown');
disp(c_hatRLS)

%-------------------------------------------------------------------------
%Remove Y from mixed audio.


% Construc echo from LMS with known delay
for m = 1:length(delay)
    for n = 1:length(DL)
        if ((n-delay(m)) > 0) &&((n-delay(m)) < length(DL))
            y_hat1(m,n)=c_hat(m)*DL(n-delay(m));
        else
            y_hat1(m,n)= 0;
        end
    end
end
Y_hat1 = sum(y_hat1);  %Sum of y_hat with different delays.
y(length(y)+1:length(Y_hat1)) = 0;
S_lms1 = s-Y_hat1';  % Mixed audio with y_hat subtracted.
E1 = ((U-S_lms1).^2); % Squared error.


%------------------------------------------------------------------------%
% Construc echo from RLS with known delay
for m = 1:length(delay)
    for n =  1:length(DL)
        if (((n-delay(m)) < length(DL)) && ((n-delay(m)) > 0))
            y_hat2(m,n)=c_hatRLS1(m)*DL(n-delay(m));
        else
            y_hat2(m,n)= 0;
        end
    end
end
Y_hat2 = sum(y_hat2);
S_rls1 = s-Y_hat2';
E2 = ((U-S_rls1).^2);


%------------------------------------------------------------------------%
% Construc echo from LMS with unknown delay
for m = 1:length(delay1)
    for n = 1:length(DL)
        if (floor(n-delay1(m)) > 0) &&((n-delay1(m)) < length(DL))
            y_hat3(m,n)=c_hatLMS(m)*DL(floor(n-delay1(m)));
        else
            y_hat3(m,n)= 0;
        end
    end
end
Y_hat3 = sum(y_hat3);
S_lms2 = s-Y_hat3';
E3 = ((U-S_lms2).^2);


%-------------------------------------------------------------------------%
% Construc echo from RLS with unknown delay
for m = 1:length(delay1)
    for n =  1:length(DL)
        if (floor(n-delay1(m)) > 0) &&((n-delay1(m)) < length(DL))
            y_hat4(m,n)=c_hatRLS(m)*DL(floor(n-delay1(m)));
        else
            y_hat4(m,n)= 0;
        end
    end
end
Y_hat4 = sum(y_hat4);
S_rls2 = s-Y_hat4';
E4 = ((U-S_rls2).^2);


t_hat = 1:length(S_lms1);  %Construct a time vector for the plots to come.


figure(4); %Plots of errors between s & y_hat.
subplot(2,2,1);
plot(t_hat,E1)
title('Known delay, LMS')
subplot(2,2,2);
plot(t_hat,E2)
title('Known delay, RLS')
subplot(2,2,3);
plot(t_hat,E3)
title('Unknown delay, LMS')
subplot(2,2,4);
plot(t_hat,E4)
title('Unknown delay, RLS')

%-------------------------------------------------------------------------%