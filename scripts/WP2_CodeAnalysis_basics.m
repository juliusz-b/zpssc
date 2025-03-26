%% WP1_CodeAnalysis.m
% This script analyzes various code families for their correlation properties
% and compares their suitability for FBG sensor interrogation systems.

% Initialize workspace
clear all;
close all;
clc;

% Ensure all folders are in the path
AddAllSubfolders;

%% Parameters
p = 8;          % Power of 2 for code length
L = 2^p-1;      % Code length
Nb = 4;         % Samples per bit
U = 8;          % Number of users/codes

%% Generate different code families
codes = struct();
codes.gold = genSpreadCodes(L, U, 'gold');
codes.kasami = genSpreadCodes(L, U, 'kasami');
codes.prbs = genSpreadCodes(L, U, 'prbs');
codes.walsh = genSpreadCodes(L, U, 'walsh');
codes.randi = genSpreadCodes(L, U, 'randi');
codes.ooc = genSpreadCodes(13, U, 'ooc',[4,2]);
codes.sidelnikov = genSpreadCodes(257, U, 'sidelnikov',[2]);
codes.golay = genSpreadCodes(L, U, 'golay');
codes.chaotic = genSpreadCodes(L, U, 'chaotic');

%%
[xc,lgs] = xcov(codes.randi(1,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.randi(2,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 60])

[xc,lgs] = xcov(codes.randi(1,:),codes.randi(2,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 60])
%%
[xc,lgs] = xcov(codes.walsh(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.walsh(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-250 250])

[xc,lgs] = xcov(codes.walsh(5,:),codes.walsh(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-100 200])

%%
[xc,lgs] = xcov(codes.prbs(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.prbs(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 60])

[xc,lgs] = xcov(codes.prbs(5,:),codes.prbs(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 60])
%%
[xc,lgs] = xcov(codes.gold(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.gold(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 60])

[xc,lgs] = xcov(codes.gold(5,:),codes.gold(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 60])

%%
[xc,lgs] = xcov(codes.kasami(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.kasami(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 60])

[xc,lgs] = xcov(codes.kasami(5,:),codes.kasami(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 60])


%%
[xc,lgs] = xcov(codes.sidelnikov(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.sidelnikov(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 60])

[xc,lgs] = xcov(codes.sidelnikov(5,:),codes.sidelnikov(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 60])


%%
[xc,lgs] = xcov(codes.ooc(1,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.ooc(2,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-10 10])
ylim([-2 5])

[xc,lgs] = xcov(codes.ooc(1,:),codes.ooc(2,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-2 5])


%%
[xc,lgs] = xcov(codes.golay(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.golay(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 260])

[xc,lgs] = xcov(codes.golay(5,:),codes.golay(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 260])
%%
[xc,lgs] = xcov(codes.chaotic(5,:));

figure('color','w');
subplot(2,1,1)
plot(lgs,xc,'LineWidth',2);
hold on;
[xc,lgs] = xcov(codes.chaotic(6,:));
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Auto kowariancja')
legend({'Kod 1','Kod 2'})
set(gca,'FontSize',15)
xlim([-30 30])
ylim([-10 30])

[xc,lgs] = xcov(codes.chaotic(5,:),codes.chaotic(6,:));
subplot(2,1,2)
plot(lgs,xc,'LineWidth',2);
xlabel('\tau')
ylabel('Kowariancja wzajemna')
set(gca,'FontSize',15)
ylim([-10 30])
