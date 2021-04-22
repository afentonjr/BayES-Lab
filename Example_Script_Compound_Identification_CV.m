% BayES' Lab Example Script
% Compound identification: this script provides an example on how to label
% compounds from an example voltammogram when provided with a library using
% cyclic voltammetry (CV).

% This file is a part of the BayES' Lab toolbox
%
% Official GitHub: """""""""""
%
% BayES' Lab: A MATLAB framework used to identify electroactive components
% from a voltammogram using a library of electrochemical-transport
% descriptors.

% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0

% Clear workspace
clear all
clc
close all

% Add to path library development directory
addpath('Compound Identification');

% Load files
PTMPTCVdatafirstelyte=load('Compound Identification\2021-03-23-PT-MPT-CV-Data-All-Trials-First-Electrolyte-Adjusted-Potential.mat');
params_indataCV=load('Compound Identification\2020-07-21-Input-Params-CV-Compound-Identify.mat');
library=load('Compound Identification\2021-04-13-Phenothiazine-Library.mat');
cv_bg=load('Compound Identification\2020-07-19-Compound-Find-CV-Background-Data-Preprocessed.mat');
cv_bgref=load('Compound Identification\2020-07-19-Compound-Find-CV-Background-Reference-Data-Preprocessed.mat');
prior=load('Compound Identification\2020-07-21-Prior-List.mat');

% Expand structure
PTMPTCVdatafirstelyte=PTMPTCVdatafirstelyte.PTMPTCVdatafirstelyte;
params_indataCV=params_indataCV.params_indataCV;
library=library.library;
cv_bg=cv_bg.cv_bg;
cv_bgref=cv_bgref.cv_bgref;
prior=prior.prior;

% Execute function
params_indataCV{7}=500; % Decrease time resolution for faster simulation
params_indataCV{9}=969.518; % Providing the correct solution resistance, in Ohms, for appropriate ohmic compensation
[Pcv,conc0cv,concfinalcv,conclbfinalcv,concubfinalcv,Estepcv,inet0cv,inetcv,dataInetcv,simtimecv]=conc_compound_find(PTMPTCVdatafirstelyte(169:end-168,1:2),params_indataCV,99,library,cv_bg,0,cv_bgref,prior,1);

% Plot probabilities as a bar graph
xlab=categorical({'PT','MPT','EPT','iPrPT','PhPT'});
xlab=reordercats(xlab,{'PT','MPT','EPT','iPrPT','PhPT'});
bar(xlab,Pcv,'LineWidth',2)
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
ylim([0 1.2])
ylabel('Posterior Probability')

% Plot the post-processed data and the best overall fit
figure()
plot(Estepcv,dataInetcv,'ob','LineWidth',2)
hold on
plot(Estepcv,inetcv,'-r','LineWidth',2)
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
xlabel('Potential (V vs. FeCp_2^{0/+})')
ylabel('Current (A)')
legend('Data','Best Fit','Location','Best')

% Plot estimated concentrations as a bar graph
xlab=xlab(Pcv>0.5); % Determine which abcissa labels should be included
figure()
bar(xlab,concfinalcv,'LineWidth',2)
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
ylabel('Estimated Concentration (mM)')