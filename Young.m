
% define physical properties
freq = 900e6; % GSM
c = 3e8;
lambda = c/freq;
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
Z0 = sqrt(mu0/epsilon0);

% define spatial and time characteristics
xdim = 200;
ydim = 200;
xsource = 75;
ysource = 100;
delta = lambda/10;
deltaT = delta/(c*sqrt(2));
Tmax = 1000; % length of the experiment

% Brain properties
eps_rel_brain = 43;
radius_brain = 10;
% Coordinates of holes
x_hole = 100;
y_hole1 = ysource + 10;
y_hole2 = ysource - 10;
z1 = x_hole + 1j*(1:y_hole2-1);
z2 = x_hole + 1j*(y_hole2+1:y_hole1-1);
z3 = x_hole + 1j*(y_hole1+1:ydim);
% build a map matrix with material relative properties
map_epsilon_rel = ones(xdim, ydim);
map_mu_rel = ones(xdim, ydim);

% build a map matrix with material absolute properties
map_epsilon = epsilon0*map_epsilon_rel;
map_mu = mu0*map_mu_rel;
map_sigma = zeros(xdim, ydim);
map_sigma_m = zeros(xdim, ydim);

% computes constant used in the FDTD equations
C_hh = (1 - (deltaT/2)* (map_sigma_m ./ map_mu)) ./ (1 + (deltaT/2)* (map_sigma_m ./ map_mu));
C_he = (deltaT ./(delta * map_mu)) ./ (1 + (deltaT/2)* (map_sigma_m ./ map_mu));
C_ee = (1 - (deltaT/2)* (map_sigma ./ map_epsilon)) ./ (1 + (deltaT/2)* (map_sigma ./ map_epsilon));
C_eh = (deltaT ./(delta * map_epsilon)) ./ (1 + (deltaT/2)* (map_sigma ./ map_epsilon));

% initialize the electric and magnetic fields component (TMz)
Hx = zeros(xdim, ydim);
Hy = zeros(xdim, ydim);
Ez = zeros(xdim, ydim);
S = zeros(xdim, ydim);
% updates the msagnetic and electric fields
% the updating corresponds to the time differentiation
Am = 1;
for t=0:1:Tmax
    % sine source
    Ez(xsource, ysource) = Am*sin(2*pi*freq*t*deltaT);
    % Maxwell equations update Hx, Hy and Ez
    for m=2:1:xdim-1
        for n=2:1:ydim-1
            Hx(m,n) = C_hh(m,n) * Hx(m,n) - C_he(m,n) * (Ez(m,n+1) - Ez(m,n));
            Hy(m,n) = C_hh(m,n) * Hy(m,n) + C_he(m,n) * (Ez(m+1,n) - Ez(m,n));
            Ez(m,n) = C_ee(m,n) * Ez(m,n) + C_eh(m,n) * ((Hy(m,n) - Hy(m-1,n)) - (Hx(m,n) - Hx(m,n-1)));
            S(m,n) = (abs(Ez(m,n))^2)/(2*Z0);
        end
    end
    % conditions of electric field on metal
    Ez(x_hole,1:1:y_hole2-1) = 0;
    Ez(x_hole,y_hole2+1:y_hole1-1) = 0;
    Ez(x_hole,y_hole1+1:ydim) = 0;
    
    %Movie type colour scaled image plot of Ez
    imagesc(1:1:xdim,(1:1:ydim)', 10*log10(abs(Ez)/Am)',[-50, 0]);colorbar; hold on;
    colormap(jet);
    plot(real(z1),imag(z1),'r'); hold on;
    plot(real(z2),imag(z2),'r'); hold on;
    plot(real(z3),imag(z3),'r');
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain at time = ',num2str(round(t*deltaT*1e+12)),' ps']); 
    xlabel('x (in cm)','FontSize',20);
    ylabel('y (in cm)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
