% Baldi Lab, 16/04/2022

% This script reads the relative permittivity data of a material and uses
% Mie theory to compute the scattering efficiency for a spherical particle
% of that material embedded in a lossless medium. The permittivity file has
% to be a tab-delimited text file with three columns: energy (in eV),
% epsilon1, epsilon2.
%
% The code then lets you manually choose at which energy you want to map
% the scattered fields, which are calculated using the equations of chapter
% 4 in Bohren and Huffman. In agreement with the book's system of
% reference, the incident plane wave propagates along the z direction, with
% the electric field polarized along the x-axis. The polar angle 'theta',
% defined as the angle between the scattered vector and the z-axis, spans
% from 0 degrees (forward scattering) to 180 degrees (back scattering), the
% azimuthal angle 'phi', defined as the angle between the x-axis and the
% projection of the scattered vector on the xy-plane, spans from 0 degrees
% to 360 degrees. The code needs the function "pin_andrea.m", corresponding
% to eq. 4.47.
%
% Finally, you have the option to plot the near-field spectrum at a
% specific location of the map. In most cases the near-field spectrum at a
% given location (Figure 3) is very different from the scattering
% efficiency of the overall particle (Figure 1).
clear all;
close all;

% MANUAL INPUT
prompt = {'Sphere radius (nm)','Refractive index of the surrounding medium (air = 1, water = 1.333)'};
dlg_title = 'Input parameters';
num_lines = 1;
def = {'',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
r=str2double(answer{1});
index=str2double(answer{2});

% DEFAULT SETTINGS
boxsize = 4*r; % Size of the 2D area in which the fields are calculated
Esteps = 200;  % Number of energy steps for the data interpolation
nmax = 10;     % Maximum order of the Bessel function
mesh = 100;    % Number of rows and columns in which the 2D box is divided (must be an even number!)
nmax = 7;      % Maximum order of the vector spherical harmonics
phi=0;         % Azimuthal angle
E0 = 1;        % Incident field intensity

% CALCULATIONS
read=dlmread(uigetfile('*','Select a dielectric function file')); % Load the permittivity file (energy in eV, eps1, eps2)
e = 1.60217646e-19; % Elementary charge in SI units
h = 6.626068e-34; % h in SI units
hbar = 1.0545718e-34; % hbar in SI units
hbareV = 6.582119514E-16; % hbar in eV*s
me = 9.10938215e-31; % Rest mass of the electron in SI units
c = 2.99792458e8; % Light speed in SI units

%Interpolates the experimental dielectric function
Emin=read(1,1);
Emax=read(size(read,1),1);
energy = (Emin:(Emax-Emin)/Esteps:Emax)';
e1_read = interp1(read(:,1),read(:,2),energy,'spline');
e2_read = interp1(read(:,1),read(:,3),energy,'spline');

%Creates the total permittivity values of the particle and the medium
n_read = (((e1_read.^2 + e2_read.^2).^(1/2) + e1_read)./2).^(1/2);
k_read = (((e1_read.^2 + e2_read.^2).^(1/2) - e1_read)./2).^(1/2);
etot_read = n_read + 1i*k_read;
m=etot_read./index; % Relative refractive index (page 100)

radius = r*1e-9; % Convertion to meters
lambda = h*c./(e*energy); % Converts energy in (eV) to wavelength in (m)
k = 2*pi*index./lambda; % Wavenumber 'k'
x=k.*radius; % Dimensionless variable (page 86)
mx=m.*x;
scaele = 0; % Initialise scattering element
extele = 0; % Initialise extinction element
for n=1:nmax
    jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
    jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
    jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
    jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
    hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
    hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
    xjnxdiff = x.*jnminx - n.*jnx;
    mxjnmxdiff = mx.*jnminmx - n.*jnmx;
    xhnxdiff = x.*hnminx - n.*hnx;
    % Calculate the scattering coefficients from eq. 4.53
    an = (m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    bn = (jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    scaele = scaele + (2*n+1).*(an.*conj(an)+bn.*conj(bn));
    extele = extele + (2*n+1).*real(an+bn);
end
Csca = 2*pi./(k.^2).*scaele; % Scattering cross section, from eq. 4.61
% Cext = 2*pi./(k.^2).*extele; % Extinction cross section, from eq. 4.62
% Cabs = Cext-Csca;            % Absorption cross section
Qsca = Csca./(pi*radius^2);  % Scattering efficiency
% Qext = Cext./(pi*radius^2);  % Extinction efficiency
% Qabs = Cabs./(pi*radius^2);  % Absorption efficiency

% Plot the scattering efficiency and choose the energy value at which the
% fields will be calculated
screen_size = get(0,'ScreenSize');
f1 = figure(1);
plot(energy,Qsca,'k','linewidth',2)
% set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
xlabel('Energy (eV)', 'FontSize', 10 );
ylabel(['Scattering efficiency of a ',num2str(radius*1E9),' nm radius sphere'],'FontSize',10);
title('Select the energy');
axis([Emin Emax 0 1.1*max(Qsca)])
[x0,y0]=ginput(1);
[dummyvalue,indexplot]=min(abs(energy-x0));
clear dummyvalue an bn x mx


% Calculate the scattered fields
E = energy(indexplot); % Energy at which we evaluate the scattered fields
lambdaplot = h*c./(e*E); % Wavelength at which we evaluate the scattered fields
kplot = 2*pi*index./lambdaplot; % Wavenumber at which we evaluate the scattered fields
x=kplot.*radius;
mx=m.*x;

% Build a 2D box in which to calculate the fields
coordx=[-boxsize/2:boxsize/mesh:boxsize/2];
coordy=[-boxsize/2:boxsize/mesh:boxsize/2];

% Initialise the fields
E2plot=zeros(size(coordx,2)/2+0.5,size(coordy,2));
% Erplot=zeros(size(coordx,2),size(coordy,2));
% Ethetaplot=zeros(size(coordx,2),size(coordy,2));
% Ephiplot=zeros(size(coordx,2),size(coordy,2));

for p=1:size(coordx,2)/2+0.5                           % I only calculate in 0<theta<180
    [num2str(round(100*p/(mesh/2+1))),'%']             % follow the progression of the calculation
    for q=1:size(coordy,2)
        distance=1e-9*sqrt(coordx(p)^2+coordy(q)^2);   % Distance from the center of the particle
        if distance>r*1e-9                             % Exclude points inside the particle
            theta=atan(abs(coordx(p))/abs(coordy(q))); % Polar angle theta
            % Adjust the polar angle theta
            if p<=(mesh+2)/2 & q<(mesh+2)/2
                theta=pi-theta;
            else if p>(mesh+2)/2 & q<(mesh+2)/2
                    theta=pi+theta;
                else if p>(mesh+2)/2 & q>=(mesh+2)/2
                        theta=2*pi-theta;
                    end
                end
            end
            costheta=cos(theta);
            
            % Initialization of variables
            Es_r = 0;
            Es_t = 0;
            Es_p = 0;
            % Hs_r = 0;
            % Hs_t = 0;
            % Hs_p = 0;
            N3_e1n_r = 0;
            N3_e1n_t = 0;
            N3_e1n_p = 0;
            M3_o1n_r = 0;
            M3_o1n_t = 0;
            M3_o1n_p = 0;
            % N3_o1n_r = 0;
            % N3_o1n_t = 0;
            % N3_o1n_p = 0;
            % M3_e1n_r = 0;
            % M3_e1n_t = 0;
            % M3_e1n_p = 0;
            
            for n=1:nmax
                
                % Prefactor, En
                En = 1i^n*E0*(2*n+1)/(n*(n+1));
                
                % Calculate the scattering coefficients at the right wavelength
                jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
                jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
                jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
                jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
                hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
                hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
                xjnxdiff = x.*jnminx - n.*jnx;
                mxjnmxdiff = mx.*jnminmx - n.*jnmx;
                xhnxdiff = x.*hnminx - n.*hnx;
                an = (m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
                bn = (jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
                
                % Hankel functions, eq. 4.13
                hnkr         = sqrt(pi/2/kplot/distance)*besselh(n+0.5,1,kplot*distance);
                hnminkr      = sqrt(pi/2/kplot/distance)*besselh(n-0.5,1,kplot*distance);
                hnpluskr     = sqrt(pi/2/kplot/distance)*besselh(n+1.5,1,kplot*distance);
                
                % Hankel function's derivative, eq. between 4.43 and 4.44
                dkrhnkr_dkr  = hnkr+kplot*distance*(n*hnminkr-(n+1)*hnpluskr)/(2*n+1);
                
                % pi_n and tau_n, eq. 4.46 and 4.47
                pin = pin_andrea(n,costheta);
                taun = (n*costheta*pin_andrea(n,costheta)-(n+1)*pin_andrea(n-1,costheta));
                
                % Components of the vector spherical harmonics, eq. 4.50
                N3_e1n_r = N3_e1n_r + cos(phi)*n*(n+1)*sin(theta)*pin*hnkr/kplot/distance;
                N3_e1n_t = N3_e1n_t + cos(phi)*taun*dkrhnkr_dkr/kplot/distance;
                N3_e1n_p = N3_e1n_p - sin(phi)*pin*dkrhnkr_dkr/kplot/distance;
                M3_o1n_r = M3_o1n_r + 0;
                M3_o1n_t = M3_o1n_t + cos(phi)*pin*hnkr;
                M3_o1n_p = M3_o1n_p - sin(phi)*taun*hnkr;
                % N3_o1n_r = N3_o1n_r + sin(phi)*n*(n+1)*sin(theta)*pin*hnkr/kplot/distance;
                % N3_o1n_t = N3_o1n_t + sin(phi)*taun*dkrhnkr_dkr/kplot/distance;
                % N3_o1n_p = N3_o1n_p + cos(phi)*pin*dkrhnkr_dkr/kplot/distance;
                % M3_e1n_r = M3_e1n_r + 0;
                % M3_e1n_t = M3_e1n_t - sin(phi)*pin*hnkr;
                % M3_e1n_p = M3_e1n_p - cos(phi)*taun*hnkr;
                
                % Components of the scattered fields, eq. 4.45
                Es_r = Es_r + En*(1i*an(indexplot)*N3_e1n_r - bn(indexplot)*M3_o1n_r);
                Es_t = Es_t + En*(1i*an(indexplot)*N3_e1n_t - bn(indexplot)*M3_o1n_t);
                Es_p = Es_p + En*(1i*an(indexplot)*N3_e1n_p - bn(indexplot)*M3_o1n_p);
                % Hs_r = Hs_r + kplot*hbareV/E*En*(1i*bn(indexplot)*N3_o1n_r + an(indexplot)*M3_e1n_r);
                % Hs_t = Hs_t + kplot*hbareV/E*En*(1i*bn(indexplot)*N3_o1n_t + an(indexplot)*M3_e1n_t);
                % Hs_p = Hs_p + kplot*hbareV/E*En*(1i*bn(indexplot)*N3_o1n_p + an(indexplot)*M3_e1n_p);
            end
            
            % Intensity of the scattered fields
            E2 = abs(Es_r).^2+abs(Es_t).^2+abs(Es_p).^2;
            % H2 = abs(Hs_r).^2+abs(Hs_t).^2+abs(Hs_p).^2;
            
            % Field maps
            E2plot(p,q) = E2;
            % H2plot(p,q) = H2;
        end
    end
end
mirror=flipud(E2plot); % Create a mirror image of the map
Emap=cat(1,E2plot,mirror([2:size(mirror,1)],:)); % Stitch the two halves together making sure to exclude the bottom row

% Plot the scattered fields
figure(2)
imagesc(log(Emap))
set(gca,'dataAspectRatio',[1 1 1])
title('Logarithm of the electric field intensity, log(|E|^2)')
xlabel('z')
ylabel('x')
button = questdlg('Would you like to plot the near-field spectrum at a specific location?')
if strcmp(button,'Yes')==1
    [qq,pp]=ginput(1);
    pp=round(pp);
    qq=round(qq);
    for p=1:size(energy)
        % Calculate the scattered fields
        E = energy(p);                  % Energy at which we evaluate the scattered fields
        lambdaplot = h*c./(e*E);        % Wavelength at which we evaluate the scattered fields
        kplot = 2*pi*index./lambdaplot; % Wavenumber at which we evaluate the scattered fields
        x=kplot.*radius;
        mx=m.*x;
        
        % Build a 2D box in which to calculate the fields
        coordx=[-boxsize/2:boxsize/mesh:boxsize/2];
        coordy=[-boxsize/2:boxsize/mesh:boxsize/2];
        
        distance=1e-9*sqrt(coordx(pp)^2+coordy(qq)^2);  % Distance from the center of the particle
        theta=atan(abs(coordx(pp))/abs(coordy(qq)));    % Polar angle theta
        % Adjust the polar angle theta
        if pp<=(mesh+2)/2 & qq<(mesh+2)/2
            theta=pi-theta;
        else if pp>(mesh+2)/2 & qq<(mesh+2)/2
                theta=pi+theta;
            else if pp>(mesh+2)/2 & qq>=(mesh+2)/2
                    theta=2*pi-theta;
                end
            end
        end
        costheta=cos(theta);
        
        % Initialization of variables
        Es_r = 0;
        Es_t = 0;
        Es_p = 0;
        N3_e1n_r = 0;
        N3_e1n_t = 0;
        N3_e1n_p = 0;
        M3_o1n_r = 0;
        M3_o1n_t = 0;
        M3_o1n_p = 0;
        
        for n=1:nmax
            
            % Prefactor, En
            En = 1i^n*E0*(2*n+1)/(n*(n+1));
            
            % Calculate the scattering coefficients at the right wavelength
            jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
            jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
            jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
            jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
            hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
            hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
            xjnxdiff = x.*jnminx - n.*jnx;
            mxjnmxdiff = mx.*jnminmx - n.*jnmx;
            xhnxdiff = x.*hnminx - n.*hnx;
            an = (m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
            bn = (jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
            
            % Hankel functions, eq. 4.13
            hnkr         = sqrt(pi/2/kplot/distance)*besselh(n+0.5,1,kplot*distance);
            hnminkr      = sqrt(pi/2/kplot/distance)*besselh(n-0.5,1,kplot*distance);
            hnpluskr     = sqrt(pi/2/kplot/distance)*besselh(n+1.5,1,kplot*distance);
            
            % Hankel function's derivative on page 94 (eq. between 4.43 and 4.44)
            dkrhnkr_dkr  = hnkr+kplot*distance*(n*hnminkr-(n+1)*hnpluskr)/(2*n+1);
            
            % pi_n and tau_n, eq. 4.46 and 4.47
            pin = pin_andrea(n,costheta);
            taun = (n*costheta*pin_andrea(n,costheta)-(n+1)*pin_andrea(n-1,costheta));
            
            % Components of the vector spherical harmonics, eq. 4.50
            N3_e1n_r = N3_e1n_r + cos(phi)*n*(n+1)*sin(theta)*pin*hnkr/kplot/distance;
            N3_e1n_t = N3_e1n_t + cos(phi)*taun*dkrhnkr_dkr/kplot/distance;
            N3_e1n_p = N3_e1n_p - sin(phi)*pin*dkrhnkr_dkr/kplot/distance;
            M3_o1n_r = M3_o1n_r + 0;
            M3_o1n_t = M3_o1n_t + cos(phi)*pin*hnkr;
            M3_o1n_p = M3_o1n_p - sin(phi)*taun*hnkr;
            
            % Components of the scattered fields, eq. 4.45
            Es_r = Es_r + En*(1i*an(p)*N3_e1n_r - bn(p)*M3_o1n_r);
            Es_t = Es_t + En*(1i*an(p)*N3_e1n_t - bn(p)*M3_o1n_t);
            Es_p = Es_p + En*(1i*an(p)*N3_e1n_p - bn(p)*M3_o1n_p);
        end
        
        % Intensity of the scattered fields
        E2 = abs(Es_r).^2+abs(Es_t).^2+abs(Es_p).^2;
        
        % Field profile
        E2spectralplot(p) = E2;
    end
    
    % Plot the scattered fields
    figure(3);
    plot(energy,E2spectralplot,'k','linewidth',2)
    title(['Electric field intensity, |E|^2, at point (',num2str(pp),',',num2str(qq),')'])
    axis([Emin Emax min(E2spectralplot) max(E2spectralplot)])
    xlabel ('Energy (eV)', 'FontSize', 10 );
end
