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
%codes.ooc = genSpreadCodes(L, U, 'ooc',[4,1]);

%%
[xc,lgs] = xcov(codes.prbs(1,:))

figure('color','w');
plot(lgs,xc);
xlabel('\tau')
ylabel('Kowariancja wzajemna')