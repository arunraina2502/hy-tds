%=============================================================================| 
% # Copyright (C) 2018 Dr.-Ing. Arun Raina (E-Mail: arunraina@icloud.com)
%
% This matlab script is part of the code used for the paper, 
% "Analysis of thermal desorption of hydrogen in metallic alloys".
% DOI: https://doi.org/10.1016/j.actamat.2017.11.011
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this code. If not, see <http://www.gnu.org/licenses/>.
%=============================================================================|  
function maps_BDH_BN
clc; clear all; close all;
format long e

bmg = 50;
NTg1 = logspace(-8,-7,10);
NTg2 = logspace(-7,-6,10);
NTg3 = logspace(-6,-5,10);
NTg4 = logspace(-5,-4,9);
NTg5 = logspace(-4,-3,8);
NTg6 = logspace(-3,-2,8);
NTg = unique([NTg1 NTg2 NTg3 NTg4 NTg5 NTg6]);
% mJi(bDH_i,bN_j) : matrix of max normalised flow rate
% mTi(bDH_i,bN_j) : matrix of peak temperature
for i = 1:bmg
    [mJi,mTi] = bDH_bN_datF(NTg(i),bmg);
    mJ(:,i) = mJi;
    mT(:,i) = mTi; i
end
dlmwrite('valuesn1.dat',mJ,'delimiter',' ','precision','%5.4d')
dlmwrite('valuesT1.dat',mT,'delimiter',' ','precision','%5.4d')

% MAIN FUNCTION =====================================================
function [mJ,mT] = bDH_bN_datF(NTa,bmg)

% define global quantities 
global Q D0 R T0 l phi NL bet NT alp DH tL0 bm

% List of input parameters
Q   = 6680;             % Lattice enthalpy [J/mol]
D0  = 2.33e-7;          % Diffusion prefactor [m2/s]
R   = 8.3144;           % Gas constant [J/K/mol]
T0  = 293;              % Temperature [K]
m   = 0;                % PDE related plane slab
l   = 5e-3;             % Specimen thickness [m]
tf  = 3.3e+4;           % Simulation final time [s]
phi = (T0*D0/l^2)*0.01; % <<<| heating rate [K/sec]
                        % [0.01 0.03 0.1]
NL  = 8.46e+28;         % NILS density [atoms/m3]      
bet = 1;                % No. of NILS per lattice atom
NT  = NL*NTa;           % Traps site density [atoms/m3]
bm  = bmg; 
bDH = linspace(10,35,bm);
DH  =-R*T0*bDH;         % Traps binding energy [J/mol]    
alp = 1;                % No. of atoms per trap site   
tL0 = 1e-6;             % Initial lat. occpancy ratio

% Space & time discretization
nx = 2e2;
nt = 3000;
xa = linspace(-l/2,l/2 ,nx);
ta = linspace(0,tf,nt);

% Normalised space and time
bQ  = Q/(R*T0); 
bphi = phi*l^2/(T0*D0);
x  = xa./l;
t  = (D0/l^2).*ta;
dx = x(2)-x(1);

% Solving PDE for \bar\theta_L at different \DeltaH 
sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
u = zeros(nt,nx,bm);
for i = 1:bm
    u(:,:,i)=sol(:,:,i);  
end

% Computing normalised flow rate at x = +l/2 at each time step
T = T0 + phi.*ta;                    % real temperature
bDL = -exp(-bQ./(1+bphi.*t))*tL0;    % normalised D_L
% 
for j = 1:bm
    for i = 1:nt
        bn(i,j) = bDL(i)*(u(i,nx,j)-u(i,nx-1,j))/dx;
    end
    findmax = false;
    s = 1;
    for i = 2:nt-1
        if T(i) > 301
        if bn(i,j)>bn(i-1,j) && bn(i,j)>bn(i+1,j) && bn(i,j)>1e-11
            mJs(j,s) = bn(i,j);
            mTs(j,s) = T(i)/T0;
            ind(j,s) = i;
            s = s+1;
            findmax = true;
%           break
        end
        else
            continue
        end
    end
    if findmax==true
        [val, ind] = max(mJs(j,:));
        mJ(j) = val;
        mT(j) = mTs(j,ind);
    elseif findmax==false
        mJ(j) = 1.0E-20;
        mT(j) = 1.0001E+00;
    end
end

% Plotting ==========================================================
% figure
% for i = 1:bm
%     plot(T./T0,bn(:,i)); hold on
%     plot(mT(i),mJ(i),'xk'); 
% end
% ylabel('J')
% xlabel('T')
% % ylim([0 3e-6])
% % xlim([290 1000])
% grid on
% set(gca,'Fontsize',13)
clear global Q D0 R T0 l phi NL bet NT alp DH tL0

% PDE coefficients ================================================== 
function [c,f,s] = pdefun(x,t,u,DuDx)
global Q D0 R T0 l phi NL bet NT alp DH tL0 bm

% Normalised quantities
bQ   = Q/(R*T0);
bphi = phi*l^2/(T0*D0);
bT   = 1 + bphi*t;
bN   = (alp*NT)/(bet*NL);
bDH  = DH./(R*T0);
K    = exp(-bDH./bT); 
bDL  = exp(-bQ/bT);
% 
for i = 1:bm
    Di(i) = 1 + (bN*K(i))/(1+tL0*K(i)*u(i))^2;
    St(i) = -(bphi*u(i))/(bT^2)*(bN*K(i)*bDH(i))/(1+tL0*K(i)*u(i))^2;
end
%
c = Di(1:bm)'; 
f = repmat(bDL,1,bm)'.*DuDx;
s = St(1:bm)';

% Initial condition ================================================= 
function u0 = icfun(x)
global bm
u0 = ones(1,bm)';

% Boundary conditions =============================================== 
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global bm
pl = ul(1:bm);
ql = zeros(1,bm)';
pr = ur(1:bm); 
qr = zeros(1,bm)';
