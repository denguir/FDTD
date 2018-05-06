% PROJECT COMMUNICATION CHANNEL
clear all;
clc;
xsource=100;
ysource=50;
xdim=200;
ydim=200;
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;
%Total no of time steps
time_tot=1000;
%stability and step
S=1/(2^0.5);
delta=1e-6;
deltat=S*delta/c;

%Epsilon matrix

Eps_mat=zeros(xdim,ydim);
x_circle=100;
y_circle=100;
Rayon=50;
for k=0:1:200
    for l=0:1:200
        if (k-100)^2+(l-100)^2 <Rayon^2
            Eps_mat(k,l)=1;
        end
    end
end

Eps=epsilon0*(ones(xdim,ydim)+Eps_mat);

% Initialization of field matrices
Ez=zeros(xdim,ydim);
%Ezx=zeros(xdim,ydim);
%Ezy=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);

%Choice of nature of source
gaussian=0;
sine=1;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=2.4e+13;
impulse=0;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step

%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                           
cst=(deltat/delta)./Eps;   
cst2=(deltat/delta)/mu0;

% Update loop begins
for n=1:1:time_tot
    %if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ez(xsource,ysource)=1;
    end
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % matrix where the wave has reached to avoid unnecessary updates.
    if n<xsource-2
        n1=xsource-n-1;
    else
        n1=1;
        %Ez(n1,n11:n21)=0;
    end
    if n<xdim-1-xsource
        n2=xsource+n;
    else
        n2=xdim-1;
        Ez(n2,n11:n21)=0;
    end
    if n<ysource-2
        n11=ysource-n-1;
    else
        n11=1;
        Ez(n1:n2,n11)=0;
    end
    if n<ydim-1-ysource
        n21=ysource+n;
    else
        n21=ydim-1;
        Ez(n1:n2,n21)=0;
    end
    
    %matrix update instead of for-loop for Ez field   
    Hx(n1:n2,n11:n21)=Hx(n1:n2,n11:n21)-cst2*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21));
    Hy(n1:n2,n11:n21)=Hy(n1:n2,n11:n21)+cst2*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21));
    Ez(n1+1:n2+1,n11+1:n21+1)=Ez(n1+1:n2+1,n11+1:n21+1)+cst(n1+1:n2+1,n11+1:n21+1).*(-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21))+cst(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1));
    
        % Source conditions
    if impulse==0
        % If unit-time step
        if gaussian==0 && sine==0
            Ez(xsource,ysource)=1;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource,ysource)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
         end
%         %if gaussian
%         if gaussian==1
%             if n<=42
%                 Ezx(xsource,ysource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/64;
%                 Ezy(xsource,ysource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/64;
%             else
%                 Ezx(xsource,ysource)=0;
%                 Ezy(xsource,ysource)=0;
%             end
%         end
    else
        %if impulse
        Ez(xsource,ysource)=0;
    end
         %Movie type colour scaled image plot of Ez
%     surf(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim)),Ez);
%     hold on;
    imagesc(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim)),Ez,[-1,1]);colorbar;
    hold on
    theta = linspace(0,2*pi,100);
    rho = 50*ones(1,100);
    z = rho.*exp(1i*theta);
    plot(real(z)+100,imag(z)+100,'b');
    hold on
    title([' Ez in a spatial domain at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
    
