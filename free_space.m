% PROJECT COMMUNICATION CHANNEL
%% all data
clear all;
clc;


frequency_ref=9e+9;
frequency=9e+8;
rho1=1.03e+3;
epsilon0=8.854197317e-12;
mu0=1.25663706e-6;
c=299792458;
sigma_e=[];
size_sigma_vec=15;
scale=4

% Setting of sigma for PML
acum=0;
for i=1:size_sigma_vec
    acum=acum+1/size_sigma_vec;
    sigma_e(i)=acum*1e-1/9;
end
sigma_m=mu0*sigma_e/epsilon0;



%Total no of time steps
time_tot=100*scale;
%stability and step
S=1/(2);
delta=(c/(10*frequency_ref));
deltat=S*delta/c;

xdim=200*scale;
ydim=200*scale;


%% Epsilon matrix

Eps_mat=ones(xdim,ydim);
mu_mat=ones(xdim,ydim);
sigma=zeros(xdim,ydim);
sigmam=zeros(xdim,ydim);




xsource=100*scale;
ysource=100*scale;


Eps=epsilon0*((Eps_mat));
Mu=mu0*(mu_mat);
% Initialization of field matrices
Ez=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);
Power=zeros(xdim,ydim);

%Multiplication factor matrices for E matrix update 

cst=((deltat/delta)./Eps)./(1+sigma*deltat./(2*Eps));   
cst2=((deltat/delta)./Mu)./(1+sigmam*deltat./(2*Mu));
cst4=(1-sigma*deltat./(2*Eps))./(1+sigma*deltat./(2*Eps));
cst5=(1-sigmam*deltat./(2*Mu))./(1+sigmam*deltat./(2*Mu));
% Update loop begins
P_profil = [];
E_profil=[];

%% plot fdtd
for n=1:1:time_tot
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
    Power(n1+1:n2+1,n11+1:n21+1)=(Ez(n1+1:n2+1,n11+1:n21+1).^2/(2*120*pi));
    
    if(n>1)
    P_profil = [P_profil Power(xsource+n,ysource)];
    E_profil= [E_profil Ez(xsource+n,ysource)];
    for k = 1:length(P_profil)
        if(Power(xsource+k,ysource) > P_profil(k))
            P_profil(k) = Power(xsource+k,ysource);
            E_profil(k)=Ez(xsource+k,ysource);
        end
    end
    end
    
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
%     hold on;
    imagesc((1:xdim)*delta,((1:ydim)*delta),10*log10(Power),[-80 0]);colorbar; %zoomed
    colormap(jet);
    hold off
    title([' Power(dB) in a spatial domain at time = ',num2str((n*sqrt(2)*deltat)),' s']); 
    xlabel('x (in m)','FontSize',20);
    ylabel('y (in m)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
%% plot
distances = (1:length(P_profil))*delta;
% figure;
% plot(distances,P_profil);
figure;
hold on
plot(distances,P_profil.*(distances),'r');

plot(distances,P_profil.*(distances.^2),'b');
legend('multiplication by r','multiplication by r^2')

figure;
plot(distances,E_profil);

% figure;
% plot(distances,1./P_profil);
% hold on
% [p,r] = polyfit(distances,1./P_profil,2);
% vec = polyval(p,distances);
% plot(distances,vec);
