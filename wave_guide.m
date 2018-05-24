% PROJECT COMMUNICATION CHANNEL
close all;
clc;
frequency_ref=9e+9;
frequency=9e+8;
rho1=1.03e+3;
rho2=1.225;
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

%Total no of time steps
time_tot=1e+8;
%stability and step
S=1/(2^0.5);
delta=(c/(10*frequency_ref));
deltat=S*delta/c;

xdim=2000;
ydim=2000;

%Required reflection co-efficient
refl_coeff=1e-6;

%Epsilon matrix
sigma_e=[];
acum=0;

size_sigma_vec=6;

for i=1:size_sigma_vec
    acum=acum+1/size_sigma_vec;
    sigma_e(i)=acum*1e-1/8;
end

sigma_m=mu0*sigma_e/epsilon0;
Eps_mat=ones(xdim,ydim);
mu_mat=ones(xdim,ydim);
sigma=ones(xdim,ydim)*3.333*1e-6;
sigmam=zeros(xdim,ydim);
Rho=rho2*ones(xdim,ydim);
x_circle=1000;
y_circle=1000;
Rayon=30;



xsource=1000;
ysource=1000;

y_wall=1:ydim;
x_wall=1050*ones(length(y_wall),1);
x_wall2=950*ones(length(y_wall),1);

sigma(x_wall,y_wall)=10000;
Eps_mat(x_wall,y_wall)=0.00001;
sigma(x_wall2,y_wall)=10000;
Eps_mat(x_wall2,y_wall)=0.00001;

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

Eps=epsilon0*((Eps_mat));
Mu=mu0*(mu_mat);
% Initialization of field matrices
Ez=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);
Power=zeros(xdim,ydim);

%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                           
cst=((deltat/delta)./Eps)./(1+sigma*deltat./(2*Eps));   
cst2=((deltat/delta)./Mu)./(1+sigmam*deltat./(2*Mu));
cst4=(1-sigma*deltat./(2*Eps))./(1+sigma*deltat./(2*Eps));
cst5=(1-sigmam*deltat./(2*Mu))./(1+sigmam*deltat./(2*Mu));
% Update loop begins
for n=1:1:time_tot
    
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
    Power(n1+1:n2+1,n11+1:n21+1)=(Ez(n1+1:n2+1,n11+1:n21+1).^2/(2*120*pi));
    sar=sigma.*(Ez.^2)./Rho;
    
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
        % If unit-time step
        %if sine
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource,ysource)=sin(((2*pi*frequency)*(n-tstart)*deltat)); 

         %Movie type colour scaled image plot of Ez
%     surf(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim)),Ez);
%     hold on;
    frac = 1/3;
    imagesc((round(frac*xdim):round((1-frac)*xdim))*delta,((round(frac*ydim):round((1-frac)*ydim))*delta),20*log10((abs(Ez((round(frac*xdim):round((1-frac)*xdim)),(round(frac*ydim):round((1-frac)*ydim)))))),[-120 0]);colorbar; %zoomed
    colormap(jet);
    hold on
    plot(y_wall*delta,x_wall*delta,'r');
    hold on
    plot(y_wall*delta,x_wall2*delta,'r');
    hold off
    title([' SAR(dB) in a spatial domain at time = ',num2str((n*sqrt(2)*deltat)),' s']); 
    xlabel('x (in m)','FontSize',20);
    ylabel('y (in m)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
    
