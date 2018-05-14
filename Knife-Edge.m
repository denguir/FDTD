% PROJECT COMMUNICATION CHANNEL
clear all;
clc;
frequency_ref=9e+9;
frequency=9e+8;
rho1=1.03e+3;
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;
sigma_e=[];
acum=0;
size_sigma_vec=6;

%rayon_head=ceil((c/frequency)*6*scale);

for i=1:size_sigma_vec
    acum=acum+1/size_sigma_vec;
    sigma_e(i)=acum*1e-1/8;
end
sigma_m=mu0*sigma_e/epsilon0;
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

Eps_mat=ones(xdim,ydim);
mu_mat=ones(xdim,ydim);
% sigma=ones(xdim,ydim)*3.333*1e-6;
sigma=zeros(xdim,ydim);
sigmam=zeros(xdim,ydim);
x_circle=1000;
y_circle=1000;
Rayon=30;

% for k=1:1:xdim
%     for l=1:1:xdim
%         if (k-x_circle)^2+(l-y_circle)^2 <(Rayon)^2
%             sigma(k,l)=1.3;
%             Eps_mat(k,l)=43;
%         end
%     end
% end
hight_source=1050;
xsource=xdim-hight_source;
ysource=450;
x_wall=xdim/2:xdim;
y_wall=1000*ones(length(x_wall),1);
sigma(x_wall,y_wall)=10000;
Eps_mat(x_wall,y_wall)=0.00001;

lambda=c/frequency;
height=abs(xsource-x_wall(1))*delta;
distance1=(y_wall(1)-ysource)*delta;
distance2=(y_wall(1)-ysource)*delta;
distance=distance1+distance2;

beta=2*pi*frequency/c;
Ei=1;
nu=height*sqrt((2/lambda)*((1/distance1) +(1/distance2)));
F=(1/2)-1i*(1/2)-fresnelc(nu)+1i*fresnels(nu);
Es=(Ei/sqrt(distance))*exp(-1i*beta*distance)*((1+1i)/2)*F;
RES=20*log10(real(Es));

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
%Ezx=zeros(xdim,ydim);
%Ezy=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);
Power=zeros(xdim,ydim);
rayons=zeros(xdim,ydim);
log10ez=zeros(xdim,ydim);


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
    SAR=sigma.*(Ez.^2)/rho1;
    Ez(x_wall,y_wall)=0;
    
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

    tstart=1;
    N_lambda=c/(frequency*delta);
    Ez(xsource,ysource)=sin(((2*pi*frequency)*(n-tstart)*deltat)); %to adapt lambda to the size of head -> frequency divided by scale
%     Ez(xsource,ysource)=1;
         %Movie type colour scaled image plot of Ez
%     surf(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim)),Ez);
%     hold on;
    imagesc((700:1300)*delta,((300:1700))*delta,20*log10(abs(Ez(700:1300,300:1700))),[-80 0]);colorbar; %zoomed
    colormap(jet);
    hold on
    plot(y_wall*delta,x_wall*delta,'r');
    hold off
    title([' SAR(dB) in a spatial domain at time = ',num2str((n*sqrt(2)*deltat)),' s']); 
    xlabel('x (in m)','FontSize',20);
    ylabel('y (in m)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
    
