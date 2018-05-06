% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "2D FDTD solution for Perfectly Matched Layer (PML) boundary condition"
% 
% Objective of the program is to solve for the Maxwell's equation for a TM 
% wave containing the xy-plane polarized magnetic field having components Hy
% and Hx and z-polarized electric field Ez. The fields are updated at every 
% timestep, in a space, where all physical parameters of free space are not
% normalized to 1 but are given real and known values. The update is done 
% using standard update equations obtained from the difference form of Maxwell's 
% curl equations with very low electric and magnetic conductivities of 4x10^(-4) 
% units incorporated in them. The field points are defined in a grid described 
% by Yee's algorithm. The H fields are defined at every half coordinate of 
% spacesteps. More precisely, the Hx part is defined at every half y-coordinate 
% and full x-coordinate and the Hy part is defined at every half x-coordinate 
% and full y-coordinate and E fields i.e the Ez part is defined at every full 
% x and full y-coordinate points.Also here, the space-step length is taken 
% as 1 micron instead of 1 unit in unitless domain assumed in previous programs. 
% Also, the time update is done using Leapfrog time-stepping. Here, H-fields
% i.e. Hx and Hy are updated every half time-step and E fields i.e Ez are 
% updated every full time-step. This is shown by two alternating matrix updates 
% spanning only a part of spatial grid where the wave, starting from source, 
% has reached at that particular time instant avoiding field updates at all 
% points in the grid which is unnecessary at that time instant. These spatial 
% updates are inside the main for-loop for time update, spanning the entire 
% time grid. Also, here, the matrices used as multiplication factors for update 
% equations are initialized before the loop starts to avoid repeated calculation 
% of the same in every loop iteration, a minor attempt at optimization. The 
% boundary condition here is Perfectly Matched Layer (PML) boundary condition 
% where the fields near the boundary are attenuated over a predetermined length 
% of boundary width before they reach the boudary to a zero value at the boundary 
% using a polynomially increasing electrical conductivity value over the boundary 
% width with maximum at the boundary and also chosing a magnetic conductivity
% value at every point in the boundary width to avoid reflection at that 
% point given in [1]. Also, here, Berenger's PML condition is used where in
% the field Ez is split into two components Ezx and Ezy and the components
% are attenuated using separate electric and magnetic conductivities in the
% two directions (sigmax and sigma_starx in x direction and sigmay and sigma_stary
% in the y direction).
%
% A source of electric field is defined at the center of the spatial domain, 
% which is a hard source, in that it does not change its value due to 
% interference from external fields i.e in other words, the source is a 
% perfect electric conductor. The form of the source can be varied using 
% the variables sine, gaussian and impulse. The source is available in four 
% standard forms- Unit-time step, Impulse, Gausian and Sinusoidal forms. 
% The color scaled plot of Ez field over the entire spatial domain is shown 
% at every time step. The simulation can be ended by closing this plot window 
% or by waiting till all the time step updates are completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing variables in memory and Matlab command screen
clear all;
clc;

% Grid Dimension in x (xdim) and y (ydim) directions
xdim=200;
ydim=200;

%Total no of time steps
time_tot=1000;

%Position of the source (center of the domain)
xsource=100;
ysource=100;

%Courant stability factor
S=1/(2^0.5);

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=1e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;

% Initialization of field matrices
Ez=zeros(xdim,ydim);
Ezx=zeros(xdim,ydim);
Ezy=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);

% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
mu=mu0*ones(xdim,ydim);

% Initializing electric conductivity matrices in x and y directions
sigmax=zeros(xdim,ydim);
sigmay=zeros(xdim,ydim);


%Perfectly matched layer boundary design
%Reference:-http://dougneubauer.com/wp-content/uploads/wdata/yee2dpml1/yee2d_c.txt
%(An adaptation of 2-D FDTD TE code of Dr. Susan Hagness)

%Boundary width of PML in all directions
bound_width=25;

%Order of polynomial on which sigma is modeled
gradingorder=6;

%Required reflection co-efficient
refl_coeff=1e-6;

%Polynomial model for sigma
sigmamax=(-log10(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width*delta);
boundfact1=((epsilon(xdim/2,bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
boundfact2=((epsilon(xdim/2,ydim-bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
boundfact3=((epsilon(bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
boundfact4=((epsilon(xdim-bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
x=0:1:bound_width;
for i=1:1:xdim
    sigmax(i,bound_width+1:-1:1)=boundfact1*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
    sigmax(i,ydim-bound_width:1:ydim)=boundfact2*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
end
for i=1:1:ydim
    sigmay(bound_width+1:-1:1,i)=boundfact3*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
    sigmay(xdim-bound_width:1:xdim,i)=boundfact4*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
end

%Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
%This is also split into x and y directions in Berenger's model
sigma_starx=(sigmax.*mu)./epsilon;
sigma_stary=(sigmay.*mu)./epsilon;

%Choice of nature of source
gaussian=0;
sine=0;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
impulse=1;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step

%Multiplication factor matrices for H matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed
G=((mu-0.5*deltat*sigma_starx)./(mu+0.5*deltat*sigma_starx)); 
H=(deltat/delta)./(mu+0.5*deltat*sigma_starx);
A=((mu-0.5*deltat*sigma_stary)./(mu+0.5*deltat*sigma_stary)); 
B=(deltat/delta)./(mu+0.5*deltat*sigma_stary);
                          
%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                          
C=((epsilon-0.5*deltat*sigmax)./(epsilon+0.5*deltat*sigmax)); 
D=(deltat/delta)./(epsilon+0.5*deltat*sigmax);   
E=((epsilon-0.5*deltat*sigmay)./(epsilon+0.5*deltat*sigmay)); 
F=(deltat/delta)./(epsilon+0.5*deltat*sigmay);

% Update loop begins
for n=1:1:time_tot
    
    %if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ezx(xsource,ysource)=0.5;
        Ezy(xsource,ysource)=0.5;
    end
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % matrix where the wave has reached to avoid unnecessary updates.
    if n<xsource-2
        n1=xsource-n-1;
    else
        n1=1;
    end
    if n<xdim-1-xsource
        n2=xsource+n;
    else
        n2=xdim-1;
    end
    if n<ysource-2
        n11=ysource-n-1;
    else
        n11=1;
    end
    if n<ydim-1-ysource
        n21=ysource+n;
    else
        n21=ydim-1;
    end
    
    %matrix update instead of for-loop for Hy and Hx fields
    Hy(n1:n2,n11:n21)=A(n1:n2,n11:n21).*Hy(n1:n2,n11:n21)+B(n1:n2,n11:n21).*(Ezx(n1+1:n2+1,n11:n21)-Ezx(n1:n2,n11:n21)+Ezy(n1+1:n2+1,n11:n21)-Ezy(n1:n2,n11:n21));
    Hx(n1:n2,n11:n21)=G(n1:n2,n11:n21).*Hx(n1:n2,n11:n21)-H(n1:n2,n11:n21).*(Ezx(n1:n2,n11+1:n21+1)-Ezx(n1:n2,n11:n21)+Ezy(n1:n2,n11+1:n21+1)-Ezy(n1:n2,n11:n21));
    
    %matrix update instead of for-loop for Ez field
    Ezx(n1+1:n2+1,n11+1:n21+1)=C(n1+1:n2+1,n11+1:n21+1).*Ezx(n1+1:n2+1,n11+1:n21+1)+D(n1+1:n2+1,n11+1:n21+1).*(-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21));
    Ezy(n1+1:n2+1,n11+1:n21+1)=E(n1+1:n2+1,n11+1:n21+1).*Ezy(n1+1:n2+1,n11+1:n21+1)+F(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1));
    
    % Source conditions
    if impulse==0
        % If unit-time step
        if gaussian==0 && sine==0
            Ezx(xsource,ysource)=0.5;
            Ezy(xsource,ysource)=0.5;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ezx(xsource,ysource)=0.5*sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
            Ezy(xsource,ysource)=0.5*sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end
        %if gaussian
        if gaussian==1
            if n<=42
                Ezx(xsource,ysource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/64;
                Ezy(xsource,ysource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/64;
            else
                Ezx(xsource,ysource)=0;
                Ezy(xsource,ysource)=0;
            end
        end
    else
        %if impulse
        Ezx(xsource,ysource)=0;
        Ezy(xsource,ysource)=0;
    end
    
    Ez=Ezx+Ezy;
    
    %Movie type colour scaled image plot of Ez
    imagesc(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim))',Ez',[-1,1]);colorbar;
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with PML boundary and at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%