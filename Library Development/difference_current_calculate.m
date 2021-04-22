function [Estep,Jnet]=difference_current_calculate(Ediff,Jdiff)
%% Header
% This function takes in the values assigned to the dimensionless flux and
% the dimensional potential from "Difference_current_extract.m" and culls
% the vector to yield the vector of baseline potentials and difference
% fluxes (to be converted to difference currents). As such, it is the
% second of two functions used to calculate the baseline potential and difference
% current. It does so by eliminating the current and potential readings
% marked for deletion by an unphysical value (e.g., 1e9).
% All values are imported from "Voltammetry_Simulator_MacroDisk_2019_11_07.m",
% and this is only used in simulating CSWVs (this processing is not needed
% to simulate conventional CVs).
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Eliminate nonsensical values from Jdiff and Ediff
Jdiff=Jdiff';
Jdiff(isnan(Jdiff)==1)=[]; % The unphysical values must match those assigned in "Difference_current_extract.m".
Ediff=Ediff';
Ediff(isnan(Ediff)==1)=[];

%% Calculate the dimensionless difference flux
% This section will determine the positive and negative dimensionless
% flux, and from there determine the dimensionless difference flux
Jdiffpos=Jdiff(1:2:end);
Jdiffneg=Jdiff(2:2:end);
Jdiffpos=Jdiffpos';
Jdiffneg=Jdiffneg';
Jnet=Jdiffpos-Jdiffneg; % This finds the dimensionless difference flux

%% Calculate the dimensional baseline potential
% This loop finds Estep, or the baseline potential
% Preallocate first
Estep=zeros(size(Ediff));
for i=1:2:length(Ediff)-1
    Estep(i)=mean([Ediff(i),Ediff(i+1)]);
    Estep(i+1)=NaN; % Assign another unrealistic value
end
    
Estep=Estep';
Estep(isnan(Estep)==1)=[]; % Cull the vector
end