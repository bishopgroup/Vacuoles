function FloryHuggins3_pHChange_1D(param)
%% FloryHuggins3_1D(param)

%% Parameters
% pH parameters
pH0 = param.pH0;
pHinf = param.pHinf;
thalf = param.thalf;
r = param.r;
[chi12,chi13,chi23] = CalcChi(pH0);

% system parameters
N = param.N;               % number of grid points (NxN) [POWER OF 2]
L = param.L;               % domain size
dt = param.dt;             % time step
A = param.A;               % numerical stability parameter
tMax = param.tMax;         % maximum time (0 < t < tmax)
Nsave = param.Nsave;       % number of configurations to be saved
plotOn = param.plotOn;     % plot results: true OR false

% thermodynamic parameters
N1 = param.N1;             % degree of polymerization of species 1 (polycation +, polymer)
N2 = param.N2;             % degree of polymerization of species 2 (polyanion -, enzyme)
N3 = param.N3;             % degree of polymerization of species 3 (solvent, water)
kappa12 = param.kappa12;   % 12 interface parameter (+-, pe)
kappa13 = param.kappa13;   % 13 interface parameter (+o, pw)
kappa23 = param.kappa23;   % 13 interface parameter (-o, ew)
phi1o = param.phi10;       % average concentration of species 1 (polycation +, polymer)
phi2o = param.phi20;       % average concentration of species 2 (polyanion -, enzyme)
phi1a = param.phi1a;       % (approximate) dilute phase concentration of species 1 (polycation +, polymer)
phi2a = param.phi2a;       % (approximate) dilute phase concentration of species 2 (polyanion -, enzyme)
phi1b = param.phi1b;       % (approximate) dense phase concentration of species 1 (polycation +, polymer)
phi2b = param.phi2b;       % (approximate) dense phase concentration of species 2 (polyanion -, enzyme)

% kinetic parameters
D11 = param.D11;           % diffusivity of species 1 (polycation +, polymer)
D22 = param.D22;           % diffusivity of species 2 (polyanion -, enzyme)
D12 = param.D12;           % cross diffusivity

%% Time points to be saved
iMax = round(tMax/dt);
tSave = linspace(0,tMax,Nsave);
iSave = round(tSave/dt);

%% Mesh
dx = L / N;
x = [0:N-1]'*dx;

%% Filenames
fileID = fopen('data/filenames.txt','w');

%% Wave numbers
l = fftshift([-N/2:N/2-1])';  % [0, 1, ..., N/2-1, -N/2, ..., -1]
kX = 2*pi*l / L;
k2 = kX.^2;
k4 = k2.^2;

%% Initial Condition 
phi1 = 0*x;
phi2 = 0*x;
gamma = 0.5*( (phi1b - phi1o)/(phi1b - phi1a) + ...
              (phi2b - phi2o)/(phi2b - phi2a) );

idx = (x < 0.5*gamma*L) | (x > L - 0.5*gamma*L);
phi1(idx) = phi1a;
phi2(idx) = phi2a;

idx = (x >= 0.5*gamma*L) & (x <= L - 0.5*gamma*L);
phi1(idx) = phi1b;
phi2(idx) = phi2b;

phi1 = smoothdata(phi1, "gaussian", 20);
phi2 = smoothdata(phi2, "gaussian", 20);
phi3 = 1 - phi1 - phi2;

PHI1 = fft(phi1);
PHI2 = fft(phi2);

SaveConfig(0,fileID,0,x,phi1,phi2,param);


%% Iterate
cntSave = 1;
for i = 2:iMax
    pH = CalcPH(dt*i,pH0,pHinf,thalf,r);
    [chi12,chi13,chi23] = CalcChi(pH);

    % Laplacian of phi
    nabla2phi1 = real(ifft(-k2.*PHI1)); 
    nabla2phi2 = real(ifft(-k2.*PHI2)); 
    nabla2phi3 = -nabla2phi1 - nabla2phi2;

    % chemical potentials
    mu1 = (1 + log(phi1)) / N1 + chi12*phi2 + chi13*phi3 + kappa12*nabla2phi2 + kappa13*nabla2phi3;
    mu2 = (1 + log(phi2)) / N2 + chi12*phi1 + chi23*phi3 + kappa12*nabla2phi1 + kappa23*nabla2phi3;
    mu3 = (1 + log(phi3)) / N3 + chi13*phi1 + chi23*phi2 + kappa13*nabla2phi1 + kappa23*nabla2phi2;
    MU1 = fft(mu1);
    MU2 = fft(mu2);
    MU3 = fft(mu3);

    % gradient of the chemical potentials
    gradXmu1 = real(ifft(1j*kX.*MU1));
    gradXmu2 = real(ifft(1j*kX.*MU2));
    gradXmu3 = real(ifft(1j*kX.*MU3));

    % diffusive flux
    j1X = D11*gFunc(phi1,N1).*(gFunc(phi3,N3) - (D12/D11)*gFunc(phi2,N2)).*(gradXmu1 - gradXmu3) ...
        + D12*gFunc(phi1,N1).*gFunc(phi2,N2).*(gradXmu2 - gradXmu3);
    j2X = D12*gFunc(phi1,N1).*gFunc(phi2,N2).*(gradXmu1 - gradXmu3) + ...
        + D22*gFunc(phi2,N2).*(gFunc(phi3,N3)- (D12/D22)*gFunc(phi1,N1)).*(gradXmu2 - gradXmu3);
    J1X = fft(j1X);
    J2X = fft(j2X);

    % Nonlinear operator
    NLO1 = 1j*kX.*J1X + A*k4.*PHI1;
    NLO2 = 1j*kX.*J2X + A*k4.*PHI2;
    
    % Subsequent Time Steps (Eq 16)
    PHI1_next = (PHI1 + NLO1*dt) ./ (1 + A*dt*k4); 
    PHI2_next = (PHI2 + NLO2*dt) ./ (1 + A*dt*k4); 
    
    % Update variables
    PHI1 = PHI1_next;
    PHI2 = PHI2_next;
    phi1 = real(ifft(PHI1));
    phi2 = real(ifft(PHI2));
    
    phi1(phi1<0) = 1e-10;
    phi2(phi2<0) = 1e-10;
    phi3 = 1 - phi1 - phi2;
    
    if i >= iSave(cntSave)
        %disp([PHI1(1),PHI2(1)])
        gamma = 0.5*( (phi1(N/2) - phi1o)/(phi1(N/2) - phi1(1)) + ...
                      (phi2(N/2) - phi2o)/(phi2(N/2) - phi2(1)) );

        % Plot results
        if plotOn == true  
            % plot(x,phi1, x,phi2); 
            semilogy(x,phi1, x,phi2); 
            % plot(x,mu1, x,mu2); 
            title(sprintf('pH = %.3f', pH),'FontSize',14);
            legend('polymer','enzyme')
            xlabel('{\itx}','FontSize',14);
            ylabel('{\it\phi}','FontSize',14);
            pause(0.1);
        end
        
        % Save Results
        SaveConfig(i,fileID,i*dt,x,phi1,phi2,param);
        cntSave = cntSave + 1;
    end
end

fclose(fileID);


function g = gFunc(phi, N)
%% g-function
g = N*phi ./ (1 + (N-1)*phi);


function SaveConfig(i,fileID,t,x,phi1,phi2,param)
%% Save results to file

% Save filename
filename = sprintf('data/data%09d.dat',i);
fprintf(fileID,'data%09d.dat\n',i);

% Generate output structure
data.t = t;
data.x = x;
data.phi1 = phi1;
data.phi2 = phi2;
data.param = param;
save(filename, '-struct','data');


function [pH] = CalcPH(t,pH0,pHinf,thalf,r)
%% return current pH

pH = pH0 + (pHinf-pH0)/(1+exp(-r*(t-thalf)));
