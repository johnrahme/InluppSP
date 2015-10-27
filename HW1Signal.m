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

%Estimation of c

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


% RLS on unknown delay
delay1 = floor([rand rand rand]*5000);
c_hatRLS = zeros(3,1);
rho = 100000;
N = 3;
P = rho*eye(N);
lambda = 0.99;
phi = zeros(3,1);
e = zeros(length(DL),1);
for n = 1:length(s)
    for m = 1:3
        if n-delay1(m)>0
            phi(m) = DL(n-delay1(m));
        else
            phi(m) = 0;
        end
    end
    %P(n) = (1/lambda)*(P(n)-(P(n)*phi*(phi')*P(n))/(lambda+((phi')*P(n)*phi)));
    P=(P-(P*phi*transpose(phi)*P)/(lambda+transpose(phi)*P*phi))/lambda;
    K = P*phi;
    
    e(n) = s(n)-c_hatRLS'*phi;
    c_hatRLS = c_hatRLS+transpose(K'*e(n));
end

disp(c_hatRLS)


% LMS on unknown delay..
delay1 = floor([rand rand rand]*5000);
c_hatLMS = zeros(3,1);
my = 0.01;
DL(length(DL)+1:length(s)) = 0;
y_hat=zeros(length(DL),1);
e=zeros(length(DL),1);
phi=zeros(3,1);
for n = 1:length(DL)
    for m = 1:3
        if n-delay1(m)>0
            phi(m)=DL(n-delay1(m));
        else
            phi(m)=0;
        end
    end
    y_hat(n)=transpose(c_hatLMS)*phi;
    e(n)=s(n)-y_hat(n);
    c_hatLMS=c_hatLMS+2*my*phi*e(n);
end


disp(c_hatLMS)








% c_hat_4 = zeros(3,1);
% data_voc(end+1:length(mixed_signl)) = 0;
% my_4 = 0.03; % Step size 0.98?
% e_4 = zeros(length(data_voc),1);
% y_hat_4 = zeros(length(data_voc),1);
% phi=zeros(3,1);
% for n = 1:length(data_voc)
%     for m=1:3
%     if n-delay(m)>0
%     phi(m)=data_voc(n-delay(m));
%     else
%      phi(m)=0;
%     end
%     end
%     y_hat_4(n) = transpose(c_hat_4)*phi;
%     e_4(n) = mixed_signl(n)-y_hat_4(n);
%     c_hat_4 = c_hat_4+2*my_4*phi*e_4(n);
% end
























