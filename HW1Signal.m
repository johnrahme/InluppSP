% Signalbehandling HW 1 %
% Vocal - near end signal
% Drum  - far  end signal
clear all
close all

[V,Fs_V] = audioread('vocal.wav');
[DL,Fs_D] = audioread('drumloop.wav');

phi = zeros(length(V),1);
theta = zeros(length(V),1);
delay = [0.05 0.1 0.3].*Fs_D; %Known delay.
c=[0.1 0.4 0.2]; %Guessed parameters. 

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
title('Echo')
subplot(2,1,2);
plot(t2,DL);
title('Far end signal')

N = sqrt(0.003)*randn(length(V),1); %Generation of random noise.
U = V+N;  %Sum of Far end signal and noise.
t_u = 1:length(U);
figure(2);
subplot(2,1,1);
plot(t_u,U)
title('Noise contaminated near end signal.')
subplot(2,1,2);
plot(t_u,V);
title('Near end signal')

%Sum of noise and echo
U(length(U)+1:length(y)) = 0;
s = y+U; %Mixed audio.
t_s = 1:length(s);
figure(3);
plot(t_s,s);
title('Audio mix')



%-------------------------------------------------------------------------
%Estimation of c usin LMS

c_hatLMS1 = zeros(3,1);
my = 0.001;  %Step size.
DL(length(DL)+1:length(s)) = 0;
y_hat=zeros(length(DL),1);
e=zeros(length(DL),1);
phi=zeros(3,1);
C_LMS_mat1 = zeros(3,length(DL));
tic; % Initiation of timer.
for n = 1:length(DL)
    for m = 1:3
        if n-delay(m)>0
            phi(m)=DL(n-delay(m));
        else
            phi(m)=0;
        end
    end
    y_hat(n)=transpose(c_hatLMS1)*phi;
    e(n)=s(n)-y_hat(n);
    c_hatLMS1=c_hatLMS1+2*my*phi*e(n);
    C_LMS_mat1(:,n) = c_hatLMS1;
end
timerLMS1 = toc; % Termination of timer.
disp('LMS known: ')
disp(c_hatLMS1);
disp(' ')

%Estimation of c using RLS
c_hatRLS1 = zeros(3,1);
%A lot of disturbance since white noise, we need a small rho.
rho = 0.01;
N = 3;
P = rho*eye(N);
lambda = 0.98;
phi1 = zeros(3,1);
e1 = zeros(length(DL),1);
C_RLS_mat1 = zeros(3,length(DL));
tic;
for n = 1:length(s)
    for m = 1:3
        if n-delay(m)>0
            phi1(m) = DL(n-delay(m));
        else
            phi1(m) = 0;
        end
    end
    K1 = P*phi1/(lambda+(phi1')*P*phi1);
    e1(n) = s(n)-c_hatRLS1'*phi1;
    c_hatRLS1 = c_hatRLS1+transpose(K1'*e1(n));    
    C_RLS_mat1(:,n) = c_hatRLS1;
end
timerRLS1 = toc;
disp('RLS known: ')
disp(c_hatRLS1)
disp(' ')


%-------------------------------------------------------------------------

% LMS on unknown delay..
N_M = 5; %upper value, when multiplied with Fs_D, in delay vector.
delay1 = (1/Fs_D:5000/Fs_D:N_M)*Fs_D; %Unknown delay vector.
c_hatLMS2 = zeros(length(delay1),1);
DL(length(DL)+1:length(s)) = 0;
y_hat=zeros(length(DL),1);
e=zeros(length(DL),1);
phi=zeros(length(delay1),1);
tic;
for n = 1:length(DL)
    for m = 1:length(delay1)
        if floor(n-delay1(m))>0
            phi(m)=DL(floor(n-delay1(m)));
        else
            phi(m)=0;
        end
    end
    y_hat(n)=transpose(c_hatLMS2)*phi;
    e(n)=s(n)-y_hat(n);
    c_hatLMS2=c_hatLMS2+2*my*phi*e(n);
end
timerLMS2 = toc;
disp('LMS unknown');
disp(c_hatLMS2)
disp(' ')

% RLS on unknown delay
c_hatRLS2 = zeros(length(delay1),1);
N = length(delay1);
P = rho*eye(N);
phi = zeros(length(delay1),1);
e = zeros(length(DL),1);
tic;
for n = 1:length(s)
    for m = 1:length(delay1)
        if floor(n-delay1(m))>0
            phi(m)=DL(floor(n-delay1(m)));
        else
            phi(m) = 0;
        end
    end
    K = P*phi/(lambda+(phi')*P*phi);
    e(n) = s(n)-c_hatRLS2'*phi;
    c_hatRLS2 = c_hatRLS2+transpose(K'*e(n));
   
end
timerRLS2 = toc;
disp('RLS unknown');
disp(c_hatRLS2)
disp(' ')
%-------------------------------------------------------------------------
%Remove Y from mixed audio.


% Construc echo from LMS with known delay
y_hat1 = zeros(length(delay1),length(DL));
for m = 1:length(delay)
    for n = 1:length(DL)
        if ((n-delay(m)) > 0) &&((n-delay(m)) < length(DL))
            y_hat1(m,n)=c_hatLMS1(m)*DL(n-delay(m));
        else
            y_hat1(m,n)= 0;
        end
    end
end
Y_hat1 = sum(y_hat1);  %Sum of y_hat with different delays.
S_lms1 = s-Y_hat1';  % Mixed audio with y_hat subtracted.
E1 = ((U-S_lms1).^2); % Squared error.


%------------------------------------------------------------------------%
% Construc echo from RLS with known delay
y_hat2 = zeros(length(delay1),length(DL));
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
y_hat3 = zeros(length(delay1),length(DL));
for m = 1:length(delay1)
    for n = 1:length(DL)
        if (floor(n-delay1(m)) > 0) &&((n-delay1(m)) < length(DL))
            y_hat3(m,n)=c_hatLMS2(m)*DL(floor(n-delay1(m)));
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
y_hat4 = zeros(length(delay1),length(DL));
for m = 1:length(delay1)
    for n =  1:length(DL)
        if (floor(n-delay1(m)) > 0) &&((n-delay1(m)) < length(DL))
            y_hat4(m,n)=c_hatRLS2(m)*DL(floor(n-delay1(m)));
        else
            y_hat4(m,n)= 0;
        end
    end
end
Y_hat4 = sum(y_hat4);
S_rls2 = s-Y_hat4';
E4 = ((U-S_rls2).^2);


t_hat = 1:length(S_lms1);  %Construct a time vector for the plots to come.


figure(4); %Plots of errors between new S and U.
subplot(2,2,1);
plot(t_hat,E1)
title('Known delay,(U-rec.signal)^2, LMS')
subplot(2,2,2);
plot(t_hat,E2)
title('Known delay,(U-rec.signal)^2, RLS')
subplot(2,2,3);
plot(t_hat,E3)
title('Unknown delay,(U-rec.signal)^2, LMS')
subplot(2,2,4);
plot(t_hat,E4)
title('Unknown delay,(U-rec.signal)^2, RLS')


%Plots of recovered signals
figure(5);
subplot(2,2,1);
plot(t_hat,S_lms1)
title('Known delay,Recovered signal, LMS')
subplot(2,2,2);
plot(t_hat,S_rls1)
title('Known delay,Recovered signal, RLS')
subplot(2,2,3);
plot(t_hat,S_lms2)
title('Unknown delay,Recovered signal, LMS')
subplot(2,2,4);
plot(t_hat,S_rls2)
title('Unknown delay,Recovered signal, RLS')

%-------------------------------------------------------------------------%
y(length(y)+1:length(Y_hat1)) = 0;
%Squared errors between y & y_hat:
E_Ylms1 = (y-Y_hat1').^2;
E_Yrls1 = (y-Y_hat2').^2;
E_Ylms2 = (y-Y_hat3').^2;
E_Yrls2 = (y-Y_hat4').^2;

figure(6); % PLots of squared errors between y & y_hat:
subplot(2,2,1);
plot(t_hat,E_Ylms1)
title('Known delay,(y-yhat)^2, LMS')
subplot(2,2,2);
plot(t_hat,E_Yrls1)
title('Known delay,(y-yhat)^2, RLS')
subplot(2,2,3);
plot(t_hat,E_Ylms2)
title('Unknown delay,(y-yhat)^2, LMS')
subplot(2,2,4);
plot(t_hat,E_Yrls2)
title('Unknown delay,(y-yhat)^2, RLS')


Y_avg_lms1 = mean((y-Y_hat1').^2);
Y_avg_rls1 = mean((y-Y_hat2').^2);
Y_avg_lms2 = mean((y-Y_hat3').^2);
Y_avg_rls2 = mean((y-Y_hat4').^2);

disp(['Time-average squared error with known delays for LMS: ', num2str(Y_avg_lms1)])
disp(['Time-average squared error with known delays for RLS: ', num2str(Y_avg_rls1)])
disp(['Time-average squared error with unknown delays for LMS: ', num2str(Y_avg_lms2)])
disp(['Time-average squared error with unknown delays for RLS: ', num2str(Y_avg_rls2)])
disp(' ')


%-------------------------------------------------------------------------%

%Plot evolution of parameters c_hat compared to real c.
t_mat=1:length(C_LMS_mat1);
c_comp = transpose(c)*ones(1, length(t_mat));
figure(7)
plot(t_mat,C_LMS_mat1(1,:),t_mat,C_LMS_mat1(2,:),t_mat,C_LMS_mat1(3,:),t_mat,c_comp)
title('Evolution of parameters using LMS')
figure(8)
plot(t_mat,C_RLS_mat1(1,:),t_mat,C_RLS_mat1(2,:),t_mat,C_RLS_mat1(3,:),t_mat,c_comp)
title('Evolution of parameters using RLS')

%-------------------------------------------------------------------------%
%Average squared errors between c & c_hat in RLS & LMS.:

%For LMS:
E_Clms = mean((c-c_hatLMS1').^2);
E_Crls = mean((c-c_hatRLS1').^2);
disp(['Average squared error in parameter using LMS: ', num2str(E_Clms)])
disp(['Average squared error in parameter using RLS: ', num2str(E_Crls)])
disp(' ')


%-------------------------------------------------------------------------%

%Time difference for parameter estimation:

disp(['Time for LMS with known delays: ', num2str(timerLMS1)])
disp(['Time for RLS with known delays: ', num2str(timerRLS1)])
disp(['Time for LMS with unknown delays: ', num2str(timerLMS2)])
disp(['Time for RLS with unknown delays: ', num2str(timerRLS2)])
disp(' ')

%-------------------------------------------------------------------------%

