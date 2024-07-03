%% prepare workspace
close all; clc; clear;

%% setup system

Rc = .14;
Rj = 1.81;
Rfj = 1.32;
Rf = .0095;
Ch = 2500;
Cc = 3800;
Cf = 1255;

C = [1 0 0
    0 0 1];
B = [1/Ch 0 0].';
A = [-1/Ch * (1/Rc + 1/Rj), 1/(Rc*Ch), 0
    1/(Rc*Cc), -1/Cc*(1/Rc + 1/Rf), 1/(Rf*Cc)
    0, 1/(Rf*Cf), -1/Cf* (1/Rfj + 1/Rf)];


T = 20;
g = ss(A,B,C,0);


%% sim
tvec = 0:T:5000;
N = length(tvec);
u = 50*[ones(1,floor(N/2)), zeros(1,floor(N/2)+1)];
Y = lsim(g,u,tvec,[0 0 0]);

figure;

plot(tvec,Y); hold on;


%% observer design

L = place(A',C',[-100, -105, -101]).'

eig(A-L*C);
