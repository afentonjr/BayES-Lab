function [Jdiff,Ediff,b]=Difference_current_extract(J,E,b,ind,maxind)
%% Header
% This function should only be used in simulating CSW voltammograms (this processing is
% not needed to simulate conventional CVs). In a CSW voltammogram, not every
% recorded potential or current is extracted - rather, for a given potential
% pulse, only a certain percentage of the full current response is averaged
% and recorded. As such, the output current needs to be processed
% accordingly. This function is the first of two to do this.
%
% In order to do this, this function takes in the previous 2 potential 
% values and decides whether these are different. If it is
% different (i.e., E(1)~=E(2)), then a new pulse has begun. As such, this function then
% “looks backwards” at the previous current values over the length
% of the pulse, averages the current values over a pre-defined fraction
% (e.g., the last 30% of the pulse), and reports that averaged
% current with the potential for that pulse. If it is not at the end of the
% pulse (i.e., E(1)=E(2)), then the function assigns an unrealistic value (e.g., 1e9) as a marker that will
% be used to cull unnecessary points during the calculation of the difference current / baseline
% potential in the function "difference_current_calculate.m". Another way
% to think about this is that the fine time mesh is being made coarse again
% (since difference currents and baseline potentials have the time
% granularity of twice the pulse duration).
% All values are imported from "Voltammetry_Simulator_MacroDisk_2019_11_07.m".
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Inputs:
%
% J -      the dimensionless flux
%
% E -      the dimensional potential (V vs. ref) of the current time point (the
%          second entry) and the previous time point (the first entry)
%
% b -      time index for the beginning of the pulse considered at the previous
%          time index
%
% ind -    the current time index
%
% maxind - the maximum index possible, to indicate the end of the
%          simulation
%
%
%% Outputs:
%
% Jdiff -   the value assigned to the dimensionless flux (either meaningful
%           or not meaningful) to calculate the difference current in "difference_current_calculate".
%
% Ediff -   the value assigned to the potential (either meaningful
%           or not meaningful) to calculate the baseline potential in "difference_current_calculate".
%
% b -       time index for the beginning of the current pulse. If the current
%           time point being investigated is the beginning of a new pulse, this value
%           will be different from the input "b".
%

%% Assign values based on the values of the two potentials
if E(2)~=E(1) % If the potential at the current time point is different from
    % that at the previous time point, then a new pulse has begun. As such,
    % this function looks backwards to calculate the difference current of
    % the previous pulse. In the settings for EC-Lab, the user can specify
    % which percentage of the current pulse is averaged. In this study, the
    % last 30% of the current response is averaged; of course, this is an
    % empirical standard and can be adjusted based on the experimental
    % parameters implemented.
    %
    % In some cases, the specified percentage may not divide neatly into
    % the number of points (e.g., when there are four points per potential pulse and the user
    % requests to average the last 30% of the current from that pulse). To address
    % this, the two indices before or after the specified percentage are extracted.
    % Average currents are calculated using these two indices (i.e.,
    % mean(current(index:end))). The average of these currents are taken in
    % turn to yield an estimate of the average current at the specified percentage.
    %
    bfrach=ceil(b+0.7.*(ind-1-b)); % Accounts for the fact that the pulse
    % averaged over the last 30% of its duration. Only counts points past 70%.
    % "h" is higher bound of points.
    %
    bfracl=floor(b+0.7.*(ind-1-b)); % Accounts for the fact that the pulse
    % averaged over the last 30% of its duration. Counts points past 70%, as
    % well as the point before 70% of the pulse. "l" is for lower bound of points
    %
    Jdiff=mean([mean(J(bfrach:ind-1)),mean(J(bfracl:ind-1))]);
    % Averaging the two definitions of how the points are extracted is
    % a first order approximation to the average of the last 30% of the
    % pulse. Naturally, the greater number of points there are, the
    % better the approximation
    %
    b=ind; % Reset counter for the next time a pulse occurs
    Ediff=E(1);
elseif ind==maxind % To account for the final point of the simulation
    bfrach=ceil(b+0.7.*(ind-1-b));
    bfracl=floor(b+0.7.*(ind-1-b));
    Jdiff=mean([mean(J(bfrach:ind-1)),mean(J(bfracl:ind-1))]);
    Ediff=E(2);
else % Assign a nonsensical value to Jdiff and Ediff, since E(2)=E(1). 
    % These points will be eliminated in "difference_current_calculate".
    Jdiff=NaN;
    Ediff=NaN;
end
end