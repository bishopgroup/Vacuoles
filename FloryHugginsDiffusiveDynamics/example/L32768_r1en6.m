%% Set Path
rootPath = pwd();
addpath(fullfile(rootPath,'../src'));

if ~exist(fullfile(rootPath,'data'), 'dir')
    mkdir(fullfile(rootPath,'data'))
end
delete 'data/*'
clear param

%% Set Parameters
% system parameters
dx = 4;
param.N = 8192;           % number of grid points (NxN) [POWER OF 2]
param.L = param.N*dx;      % domain size
param.dt = 5;              % time step
param.A = 40;              % numerical stability parameter
param.Nsave = 1000;        % number of configurations to be saved
param.plotOn = false;      % plot results: true OR false

% pH parameters
param.pH0 = 10;
param.pHinf = 8;
param.tMax = 6 * 0.01 * param.L^2;  
param.thalf = 2 * 0.01 * param.L^2; 
param.r = 1e-6;
[chi12,chi13,chi23] = CalcChi(param.pH0);
Delta = (chi13 - chi23) / (chi13 + chi23 - chi12 - 1);
alpha = 1/(Delta + sqrt(1 + Delta^2));

% thermodynamic parameters
param.N1 = 10704.8;        % degree of polymerization of species 1 (polycation +, polymer)
param.N2 = 2709.33;        % degree of polymerization of species 2 (polyanion -, enzyme)
param.N3 = 1;              % degree of polymerization of species 3 (solvent, water)
param.kappa12 = 1;         % 12 interface parameter (+-, pe)
param.kappa13 = 0.5;       % 13 interface parameter (+o, pw)
param.kappa23 = 0.5;       % 23 interface parameter (-o, ew)
param.phi10 = 0.003;                 % average concentration of species 1 (polycation +, polymer)
param.phi20 = alpha * param.phi10;   % average concentration of species 2 (polyanion -, enzyme)
param.phi1a = 0.000280451;      % (approximate) dilute phase concentration of species 1 (polycation +, polymer)
param.phi2a = 0.000335796;      % (approximate) dilute phase concentration of species 2 (polyanion -, enzyme)
param.phi1b = 0.0286586;        % (approximate) dense phase concentration of species 1 (polycation +, polymer)
param.phi2b = 0.0316921;        % (approximate) dense phase concentration of species 2 (polyanion -, enzyme)

% kinetic parameters
param.D11 = 1;             % diffusivity of species 1 (polycation +, polymer)
param.D22 = 1;             % diffusivity of species 2 (polyanion -, enzyme)
param.D12 = 0;             % cross diffusivity

%% Run Flory-Huggins
tic
FloryHuggins3_pHChange_1D_fast(param);
toc

%% Clean path
rmpath(fullfile(rootPath,'../src'));
