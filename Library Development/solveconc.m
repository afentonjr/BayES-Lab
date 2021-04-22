function conc=solveconc(X,params,theta,cinit,h)
%% Header
% This function is the solver for the concentration profile of a
% voltammogram at a particular point in time during a voltammetric
% simulation. The method used here is a finite difference stencil applied to
% Fick's Second Law in one spatial dimension. This solution method is
% very similar to that employed in Understanding Voltammetry: Simulation of
% Electrode Processes by Compton, Laborda, and Ward (Imperial College Press,
% 2014); specifically Chapters 3.2.2-3.4.1 and 4.1-4.3.
%
% This function takes in the dimensionless spatial mesh, the dimensionless overpotential, and
% the initial concentration guess (which, in a voltammametric simulator,
% the vector of concentrations at previous time points), and necessary parameters.
% The function then outputs the vector of concentrations at the given
% dimensionless overpotential by implicitly solving the linear algebra problem
% Ax=b, where A is a matrix of constants, x is the vector of unknown concentrations,
% and b is the vector of initial concentrations.
% All parameter values are imported from the parent function "Voltammetry_Simulator_MacroDisk_2019_11_07.m".
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Inputs:
%
% X -      the spatial mesh, dimensionless.
%
% params - imported parameters. The order of the parameters is:
%          (1) alpha, the transfer coefficient (unitless)
%          (2) Ko, the dimensionless heterogeneous rate constant
%          (3) dtau, the dimensionless time step
%          (4) d_B, the ratio of the two diffusion coefficients D_B/D_A
%          (5) nxpoints, the number of spatial points
%          (6) c_A_inf, the dimensionless bulk concentration of A
%          (7) c_B_inf, the dimensionless bulk concentration of B
%          (8) kinetics, the electron transfer mechanism being considered (1 = reversible, other number = quasireversible).
% 
% theta -  the dimensionless overpotential
%
% cinit -  the initial dimensionless concentration of both A and B. In a voltammetric
%          simulator, it is the vector of concentrations from the previous time point.
%
% h -      the initial step size (i.e., X(2)-X(1)). Dimensionless
%
%
%% Output:
%
% conc -   the dimensionless concentration vector. The ordering of the
%          entries with respect to the spatial coordinate is as follows:
%          [c_A(max_distance), c_A(max_distance - 1),...,c_A(2),c_A(1),c_B(1),c_B(2),...,
%          c_B(max_distance-1),c_B(max_distance)].

%% Extract parameters
% Params
alpha = params(1);
Ko = params(2);
dtau = params(3);
d_B = params(4);
nxpoints = params(5);
c_A_inf = params(6);
c_B_inf = params(7);
kinetics = params(8);

%% Set up surface boundary points
xdivide=nxpoints-2;
delt=cinit; % In this case, "delt" is the right-hand side of the equation Ax=b,
% where "A" is a matrix of coefficients, "x" is the vector of
% concentrations being solved for, and "b" is as described above.
% Electrode boundary conditions
if kinetics==1 % Reversible
    delt(end/2)=0;
    delt(end/2+1)=0;
    alphaA(nxpoints)=0;
    gammaA(nxpoints)=0;
    alphaB(1)=0;
    gammaB(1)=0;
    betaA(nxpoints)=1;
    betaB(1)=1;
else % Quasireversible
    ftheta=exp(-alpha.*theta);
    alphaA(nxpoints)=-h*ftheta*Ko; % Alpha scale (i.e., of species initially present)
    betaA(nxpoints)=1+h*ftheta*Ko*exp(theta); % Beta scale
    gammaA(nxpoints)=-1;
    delt(end/2)=0;
    alphaB(1)=-(h/d_B)*exp(theta)*ftheta*Ko; % Alpha not scale (i.e., of species not initially present)
    betaB(1)=1+(h/d_B)*ftheta*Ko;
    gammaB(1)=-1;
    delt(end/2+1)=0;
end

%% Set up interior points
for i=2:xdivide+1
    alphaA(nxpoints+1-i)=-2*dtau/((X(i+1)-X(i-1))*(X(i)-X(i-1)));
    gammaA(nxpoints+1-i)=-2*dtau/((X(i+1)-X(i-1))*(X(i+1)-X(i)));
    betaA(nxpoints+1-i)=-alphaA(nxpoints+1-i)-gammaA(nxpoints+1-i)+1;
    alphaB(i)=-2*dtau*d_B/((X(i+1)-X(i-1))*(X(i)-X(i-1)));
    gammaB(i)=-2*dtau*d_B/((X(i+1)-X(i-1))*(X(i+1)-X(i)));
    betaB(i)=-alphaB(i)-gammaB(i)+1;
end

%% Set up infinite BCs
alphaA(1)=0;
betaA(1)=1;
delt(1)=c_A_inf;
alphaB(nxpoints)=0;
betaB(nxpoints)=1;
delt(2*nxpoints)=c_B_inf;
gammaA=gammaA(2:end);

alpha=[gammaA,alphaB];
beta=[betaA,betaB];
gamma=[alphaA,gammaB];

%% Solve for the concentration vector using linear algebra
if kinetics==1 % Nernstian kinetics
    C=diag(alpha,-1)+diag(beta)+diag(gamma,1);
    % Non-first-subdiagonal entries
    C(nxpoints,nxpoints-1:nxpoints+2)=[(1+d_B*exp(theta))^(-1),-1,0,d_B/(1+d_B*exp(theta))];
    C(nxpoints+1,nxpoints-1:nxpoints+2)=[0,exp(theta),-1,0];
    delt(nxpoints)=0;
    delt(nxpoints+1)=0;
    C=sparse(C);
    conc=C\delt;
else % Butler Volmer kinetics. Backsolving is used here, rather than the backslash operator,
    % because less memory is needed to store the necessary information
    % (decreasing simulation time). It is important to note that
    % this method (i.e., the Thomas Algorithm) only works for tridiagonal
    % matrices, which the Nernstian kinetics case does not satisfy
    gammamod(1)=gamma(1)/beta(1);
    deltmod(1)=delt(1)/beta(1);
    for i=2:2*nxpoints-1
        gammamod(i)=gamma(i)/(beta(i)-gammamod(i-1)*alpha(i-1));
        deltmod(i)=(delt(i)-deltmod(i-1)*alpha(i-1))/(beta(i)-gammamod(i-1)*alpha(i-1));
    end
    deltmod(2*nxpoints)=(delt(2*nxpoints)-deltmod(2*nxpoints-1)*alpha(2*nxpoints-1))/...
        (beta(2*nxpoints)-gammamod(2*nxpoints-1)*alpha(2*nxpoints-1));
    conc(2*nxpoints)=deltmod(2*nxpoints);
    for i=2*nxpoints-1:-1:1
        conc(i)=deltmod(i)-gammamod(i)*conc(i+1);
    end
end
end