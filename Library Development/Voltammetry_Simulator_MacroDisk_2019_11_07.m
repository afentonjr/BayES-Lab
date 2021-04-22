function [conc,potential,I,t,phiLSE]=Voltammetry_Simulator_MacroDisk_2019_11_07(kinetics,direction,type,dataV,params)
%% Header
% Simulates a CV or CSWV of a redox couple: A <=> B^+ + e- on a
% macroelectrode disk (e.g., 1.5 mm radius glassy carbon working
% electrode), given knowledge of the kinetic mechanism (reversible or
% quasireversible), the direction of the initial sweep, the type of
% voltammogram (CV or CSWV), the input potential data, and necessary
% paramters.
% Assumes full support (infinite supporting electrolyte)
% A is ALWAYS reduced species, B is ALWAYS oxidized species
%
% Note that this simulator works if there is only one species (either
% oxidized or reduced) in the bulk. If there is ambiguity as to whether
% there is only oxidized or reduced species in the bulk, then a low or
% high-potential hold can be performed to make the bulk concentration
% entirely reduced or oxidized, respectively.
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Inputs:
% kinetics -  a counter for whether a reversible (kinetics = 1) or quasireversible
%             (kinetics = any other number) voltammogram is simulated
%
% direction - an indicator for whether the initial sweep is oxidative or
%             reductive. If direction = 1, the initial sweep is oxidative; otherwise,
%             it is reductive.
%
% type -      simulate either a conventional CV or CSWV. If type = 1, then a
%             conventional CV is simulated. Otherwise, a CSWV is simulated.
%
% dataV -     a vector of potentials in the waveform of the desired voltammetric
%             technique (i.e., linear waveform for CV, square waves for CSWV), fed from
%             the program that is using this voltammetric simulator. In units of V vs. Ref
%
% params -    the relevant parameters for the simulation. The order of the
%             parameters is
%             (1) the bulk concentration in mol/L
%             (2) the diffusion coefficient of the reduced species in cm^2/s
%             (3) the diffusion coefficient of the oxidized coefficient in cm^2/s
%             (4) the radius of the working electrode in mm
%             (5) the step height of the CSWV pulse in V
%             (6) the pulse height of the CSWV in V
%             (7) the pulse duration of the CSWV in s
%             (8) the temperature of the solution in K
%             (9) the formal potential in V vs. ref
%             (10) the transfer coefficient (unitless, if applicable)
%             (11) the heterogeneous rate constant in cm/s (if applicable)
%             (12) the potential step of the simulation in V
%             (13) the scan rate in V/s
%
%             If a reversible electron transfer is being simulated, arbitrary
%             values of the transfer coefficient and the heterogeneous rate constant
%             may be used.
%
%
%% Outputs:
%
% conc -      the concentration matrix, in mol/m^3. Columns are the
%             concentration as a function of space. Rows are the concentration at a
%             single spatial point over time.
%
% potential - a vector of potentials, one for each time point. In V vs. ref.
%
% I -         a vector of the output current (either the raw current or difference
%             current), one for each time point. In A.
%
% t -         a vector of the simulation time, in s.
%
% phiLSE -    the concentration-normalized current, used for fitting the
%             library to simulated voltammograms. In units of A*m^3/mol.

%% Extract parameters
% A is ALWAYS the reduced species, B is ALWAYS the oxidized species. The
% scales can vary though, depending on the direction of the sweep
C_inf=params(1);
D_A=params(2); % ALWAYS Dr
D_B=params(3); % ALWAYS Do
% The below allows for the flexibility to simulate CVs in either direction
if direction==1 % From negative to positive potentials initially (i.e., only reduced species present initially)
    C_A_inf = C_inf; % initial concentration of A, the reduced species (mol/L)
    C_B_inf = 0; % initial concentration of B, the oxidized species (mol/L)
    C_scale=(C_A_inf*1000); % Cscale is Cr, convert to mol/m^3
    [~,indmax]=max(dataV);
    if type==1 % CV
        Eo=dataV(1); % Starting potential
        Elam=dataV(indmax); % The switching potential is the maximum potential
    else % CSWV
        Eo=mean([dataV(1),dataV(2)]); % Starting baseline potential is the mean of the first two potential points
        Elam=mean([dataV(indmax),dataV(indmax+1)]); % The counter after the maximum potential is still in the forward sweep
    end
else % From positive to negative potentials initially (i.e., only oxidized species present initially)
    C_B_inf = C_inf; % initial concentration of B, the oxidized species (mol/L)
    C_A_inf = 0; % initial concentration of A, the reduced species (mol/L)
    C_scale=(C_B_inf*1000); % Cscale is Co
    [~,indmin]=min(dataV);
    if type==1
        Eo=dataV(1);
        Elam=dataV(indmin);
    else
        Eo=mean([dataV(1),dataV(2)]);
        Elam=mean([dataV(indmin),dataV(indmin+1)]);
    end
end
re=params(4); % electrode radius (mm)
dE=params(5); % CSWV, V per step
Esw=params(6); % CSWV, pulse amplitude (V)
tpulse=params(7); % CSWV, pulse duration (s)
T=params(8); % K, solution temperature
Ef=params(9); % (V vs. ref, formal redox potential of the couple)
alpha=params(10); % (unitless, transfer coefficient)
ko=params(11); % (cm/s, heterogeneous rate constant)
dnu=params(12); % Discretization of applied potential (deltaV per simulated point)
if type==1 % If a conventional CV (cyclic linear sweeep voltammetry) is simulated
    nu=params(13); % Directly use the value of the scan rate, in mV/s
else % IF a CSWV is being simulated
    nu=(dE/2/tpulse)*1000; % Adjusted scan rate. Converted back to mV/s
end

%% Convert to SI units
C_A_inf=C_A_inf.*1000; % mol/m^3
C_B_inf=C_B_inf.*1000; % mol/m^3
nu=nu./1000; % V/s
re=re./1e3; % m
D_A=D_A./1e4; % cm^2/s to m^2/s
D_B=D_B./1e4; % cm^2/s to m^2/s
ko=ko./100; % cm/s to m/s

%% Unchanging parameters
F=96485; % C/mol
Rgas=8.314; % J/mol/K

%% Determine time mesh based on potential
% The time steps between each of the potential values are of a certain granularity
% (e.g., in CSWV, it is that of the pulse duration). In the case of CSWV,
% if the pulse duration is one second, then the time step between the
% reported potentials is also one second. This time step may be be too
% large to accurately simulate voltammograms, so the time step size usually
% must be decreased. This increases the number of data points, and hence
% simulation time, but also increases the accuracy of the simulator. The argument
% is similar for a conventional CV, even if the original time mesh is finer. For both techniques,
% the current and potential are later transformed and output with a more coarse discretization
% (i.e., the time steps are made coarse again) to match the granularity of the
% input data for consistent comparison.
dt=dnu./abs(nu); % Desired time step (s), based on the desired potential discretization
sub=2*abs(Elam-Eo)/dnu; % The number of substeps
sub=ceil(sub); % Make it an integer
if type==1 % The input time metric differs depending on whether a CV or a CSWV is simulated
    [t,E]=E_find_interpolate(dataV,sub,dt*sub/length(dataV)); % Expand the time mesh to have "sub" number of substeps (or time points).
    % The time step for CVs needs an additional factor (sub/length(dataV))
    % for the math to work out correctly.
else
    [t,E]=E_find_interpolate(dataV,sub,tpulse);
end
tmax=max(t);
% Non-dimensionalize t to be tau
tau=D_A.*t./re.^2;
dtau=D_A.*dt./re.^2;
% Non-dimensionalize E to be theta
theta=(F./Rgas./T).*(E-Ef);
%% Initialize spatial mesh in x-direction (Cartesian coordinates)
m=1; % initialize counter
h=1e-5.*re; % Distance between first two points
omega=1.1; % Exponential factor for x
limdist=6.*sqrt(max([D_A,D_B]).*tmax); % Distance of simulation, tmax is the maximum time
X(1)=0; % Initialize distance variable
while X(m)<=limdist % Determine distance
    X(m+1)=h.*omega.^(m-1)+X(m); % Determine the next spatial point to be studied (exponentially expanding mesh)
    m=m+1; % Update counter
end
nxpoints=m; % Length of the radius spatial vector
% Non-dimensionalize x
h0=h./re;
X=X./re;
%% Non-dimensionalize other variables
% Concentrations
c_A_inf=C_A_inf./C_scale;
c_B_inf=C_B_inf./C_scale;

% Diffusion coefficients
d_B=D_B./D_A;

% Heterogeneous rate constant
K0=ko.*re./D_A;
%% Initial concentrations
cinit=[flip(ones(nxpoints,1).*c_A_inf);ones(nxpoints,1).*c_B_inf];

%% Create parameter vector to be passed into the solver function
params(1) = alpha;
params(2) = K0;
params(3) = dtau;
params(4) = d_B;
params(5) = nxpoints;
params(6) = c_A_inf;
params(7) = c_B_inf;
params(8) = kinetics;
%% Solve for concentrations over time
b=1; % Index that marks the beginning of the current potential pulse
for i=1:length(tau)
    c=solveconc(X,params,theta(i),cinit,h0);
    conc(:,i)=c; % Export the concentration at a given time point to the output concentration matrix
    cinit=c; % Initialize the concentration for the next time step
    if direction==1 % The current is defined with respect to that of species A (the reduced species).
        % The indices for the reduces species depend on the direction of
        % the initial sweep (i.e., whether it is oxidative or reductive).
        % As a result, the extracted indices from the concentration vector
        % differ depending on the direction. In the case of an initially
        % oxidative sweep, the indices below are extracted.
        c2=c(end/2-1);
        c1=c(end/2);
    else % If the initial sweep direction is reductive, the below indices are extracted.
        c2=c(end/2+2);
        c1=c(end/2+1);
    end
    J(i)=(c2-c1)./h0; % Find the flux, which is uniform on the surface
    if i>1 % Creates a vector of difference currents and baseline potentials for further processing
        [Jdiff(i-1),Ediff(i-1),b]=Difference_current_extract(J,E(i-1:i),b,i,length(tau));
    end
end
if type~=1 % Culls and calculates the difference current as a function of baseline potential
    [E,J]=difference_current_calculate(Ediff,Jdiff);
else % Compresses the potential and flux to the original granularity of the CV
    Ediff(isnan(Ediff)==1)=[];
    E=Ediff;
    Jdiff(isnan(Jdiff)==1)=[];
    J=Jdiff;
end

%% Define currents over all potentials
phiLSE=sign(Elam-Eo)*J.*pi.*re.*F.*D_A; % For least-squares fitting when the library is compared to experimental data
I=phiLSE*C_scale; % Redimensionalize current
potential=E; % Potential vector
end