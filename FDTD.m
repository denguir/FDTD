
% define physical properties
freq = 900e6; % GSM
freq_ref = 10*freq;
c = 3e8;
lambda = c/freq;
lambda_ref = c/freq_ref;
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
Z0 = sqrt(mu0/epsilon0);

% define spatial and time characteristics
xdim = 400;
ydim = 400;
xsource = 200;
ysource = 200;
delta = lambda_ref/10;
deltaT = delta/(c*sqrt(2));
Tmax = 1000; % length of the experiment

% Brain properties
eps_rel_brain = 43;
radius_brain = 30; % 30pixel = 20cm (scale factor due to lambda_ref)
% Coordinates of brain
xbrain = xsource + radius_brain;
ybrain = ysource;
theta = linspace(0,2*pi,100);
rho = radius_brain*ones(1,100);
z = rho.*exp(1i*theta);

% build a map matrix with material relative properties
map_epsilon_rel = ones(xdim, ydim);
map_mu_rel = ones(xdim, ydim);
for x=1:xdim
    for y=1:ydim
        if (x-xbrain)^2 + (y-ybrain)^2 < radius_brain^2
            map_epsilon_rel(x,y) = eps_rel_brain;
        end
    end
end

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
    
    %Movie type colour scaled image plot of Ez
    imagesc(1:1:xdim,(1:1:ydim)', 20*log10(abs(Ez)/Am)',[-100, 0]);colorbar; hold on;
    colormap(jet);
    plot(real(z)+xbrain,imag(z)+ybrain,'r');
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain at time = ',num2str(round(t*deltaT*1e+12)),' ps']); 
    xlabel('x (in delta_unity)','FontSize',20);
    ylabel('y (in delta_unity)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
