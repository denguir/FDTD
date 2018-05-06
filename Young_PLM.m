% PROJECT COMMUNICATION CHANNEL
clear all;
clc;

frequency=7e+9;
rho1=1.03e+3;
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;
scale=10^(-floor(log10(c/frequency)));
sigma_e=[];
acum=0;
size_sigma_vec=6;

for i=1:size_sigma_vec
    acum=acum+1/size_sigma_vec;
    sigma_e(i)=acum*1e-1/8;
end
sigma_m=mu0*sigma_e/epsilon0;
%Total no of time steps
time_tot=1e+8;
%stability and step
S=1/(2^0.5);
delta=c/(10*frequency);
deltat=S*delta/c;
xsource=100;
ysource=100;
xdim=200;
ydim=200;

%Required reflection co-efficient
refl_coeff=1e-6;

%Epsilon matrix
Eps_mat=zeros(xdim,ydim);
mu_mat=zeros(xdim,ydim);
sigma=ones(xdim,ydim)*3.333*1e-6;
sigmam=zeros(xdim,ydim);

% Coordinates of holes
y_hole = ysource+25;
x_hole1 = xsource + 10;
x_hole2 = xsource - 10;
z1 = y_hole + 1j*(1:x_hole2-1);
z2 = y_hole + 1j*(x_hole2+1:x_hole1-1);
z3 = y_hole + 1j*(x_hole1+1:xdim);


for i=1:length(sigma_e)
    sigma(i,:)=sigma_e(i);
    sigma(xdim-size_sigma_vec-1+i,:)=sigma_e(i);
    sigma(:,i)=sigma_e(i);
    sigma(:,ydim-size_sigma_vec-1+i)=sigma_e(i);
end
for i=1:length(sigma_m)
    sigmam(i,:)=sigma_m(i);
    sigmam(xdim-size_sigma_vec-1+i:xdim-1,:)=sigma_m(i);
    sigmam(:,i)=sigma_m(i);
    sigmam(:,ydim-size_sigma_vec-1+i)=sigma_m(i);
end
Eps=epsilon0*(ones(xdim,ydim)+(Eps_mat));
Mu=mu0*(ones(xdim,ydim)+mu_mat);
% Initialization of field matrices
Ez=zeros(xdim,ydim);
%Ezx=zeros(xdim,ydim);
%Ezy=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);
Power=zeros(xdim,ydim);
log10ez=zeros(xdim,ydim);
   
%Choice of nature of source
sine=1;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
impulse=0;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step

%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                           
cst=((deltat/delta)./Eps)./(1+sigma*deltat./(2*Eps));   
cst2=((deltat/delta)./Mu)./(1+sigmam*deltat./(2*Mu));
cst3=((c/frequency)^2)/((2*120*pi*4*pi));
cst4=(1-sigma*deltat./(2*Eps))./(1+sigma*deltat./(2*Eps));
cst5=(1-sigmam*deltat./(2*Mu))./(1+sigmam*deltat./(2*Mu));
% Update loop begins

for n=1:1:time_tot
        %if source is impulse or unit-time step 
    if  sine==0 && n==1
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
    
    %matrix update instead of for-loop for Ez field  

    
    Hx(n1:n2,n11:n21)=cst5(n1:n2,n11:n21).*Hx(n1:n2,n11:n21)-cst2(n1:n2,n11:n21).*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21));
    Hy(n1:n2,n11:n21)=cst5(n1:n2,n11:n21).*Hy(n1:n2,n11:n21)+cst2(n1:n2,n11:n21).*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21));
    
    Ez(n1+1:n2+1,n11+1:n21+1)=cst4(n1+1:n2+1,n11+1:n21+1).*Ez(n1+1:n2+1,n11+1:n21+1)+cst(n1+1:n2+1,n11+1:n21+1).*(-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21))+cst(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1));
    Power(n1+1:n2+1,n11+1:n21+1)=cst3*(Ez(n1+1:n2+1,n11+1:n21+1).^2);
    log10ez=mag2db(Ez);
    SAR=sigma.*(Ez.^2)/rho1;
    
    % conditions of electric field on metal
    Ez(1:1:x_hole2-1,y_hole) = 0;
    Ez(x_hole2+1:x_hole1-1,y_hole) = 0;
    Ez(x_hole1+1:xdim, y_hole) = 0;
    
    if n11==1
        Ez(:,n11)=0;
    end
    
    if n21==ydim-1
        Ez(:,n21)=0;
    end
    
    if n1==1
        Ez(n1,:)=0;
    end
    
    if n2==xdim-1
        Ez(n2,:)=0;
    end
    % Source conditions
    if impulse==0
        % If unit-time step
        if sine==0
            Ez(xsource,ysource)=1;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource,ysource)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end

    else
        %if impulse
        Ez(xsource,ysource)=0;
    end
         %Movie type colour scaled image plot of Ez
%     surf(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim)),Ez);
%     hold on;
    imagesc((1:xdim),((1:ydim)'),20*log10(abs(Ez)),[-100 0]);colorbar;
    plot(real(z1),imag(z1),'r'); hold on;
    plot(real(z2),imag(z2),'r'); hold on;
    plot(real(z3),imag(z3),'r');
    title([' Power (dBm) in a spatial domain at time = ',num2str(round(n*deltat*1e9)),' ns']); 
    xlabel('x (in mm)','FontSize',20);
    ylabel('y (in mm)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
    
