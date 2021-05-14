

clear all
close all
%
M=10;          % sensor number
T=20;          % snapshot number
SNR=2;        % SNR
tt=1:T;
L=length(tt);
error=0;error1=0;
test=200; % test number
timetol=0;iter=0;timetol3=0;timetol1=0;iter1=0;
for i=1:test  % independent simulation
grid_interval=2;  % Azimuth grid_interval
Azimuth_grid=-90:grid_interval:90;  % all hypothetical DOA angles 
%% random DOAs£º
u=unifrnd(-grid_interval,grid_interval);
AL=[ -7+u 3+u 10+u ]; 
%% fixed DOAs
% AL=[ -7.05 3.57 10.27 ];  
K=size(AL,2);  % signal number     
%% ULA array manifold
A=[];
 for m=1:M
    for i=1:K 
    A(m,i)=exp(1j*pi*(m-1)*sin(AL(1,i)*pi/180));
    end
 end
 Aa=[];
 %% the ULA array completed manifold matrix:
    for m=1:M    
    for ang= 1:length(Azimuth_grid)
    Az(ang)=Azimuth_grid(ang)*pi/180;
    Aa(m,ang)=exp(1j*pi*(m-1)*sin(Az(ang)));
    end
    end
    

%% 1.Uncorrelated stochastic signal :
S=(randn(K,L)+1j*randn(K,L));  % complex signals          
Vj=diag(sqrt(  10^(SNR/10)./diag(1/L*(S*S') ) ) ); 
S=Vj*S;     
noise=sqrt(1/2)*(randn(M,L)+1j*randn(M,L));  % complex Gaussian noise
Y=A*S+noise;

%% 2 Correlated  stochastic signals :
%  amp = [1 1]';
% x1 = (randn(1,L) + 1j *randn(1,L))/sqrt(2);
% y1 = (randn(1,L) + 1j *randn(1,L))/sqrt(2);
% cor=1;
% x2 = cor*x1 + sqrt(1-cor^2)*y1;
% S = diag(amp) * [x1; x2];
% Vj=diag(sqrt(  10^(snr/10)./diag(1/L*(S*S') ) ) ); 
% S=Vj*S;     %
% noise=sqrt(1/2)*(randn(M,L)+1j*randn(M,L));  % 
% Y=A*S+noise;

%%
tstart1 = tic;
[Source_power_re,Azimuth,source_id_,j1]=CGDP_SBL(Y,Azimuth_grid,grid_interval,Aa,K);   % K is known
% [Source_power_re,Azimuth,source_id_,j1]=CGDP_FSBL(Y,Azimuth_grid,grid_interval,Aa,M-1);  % if K is unknown
DOA_estimate=Azimuth(source_id_);
time1=toc(tstart1);
timetol1=timetol1+time1  %total time
iter1=iter1+j1 % iteration number
error_sum1=0;
 if length( DOA_estimate)>=K
   for i=1:length(AL)
    error_sum1=error_sum1+((DOA_estimate(i)-AL(i)))^2;
   end
 else 
    for i=1:length(DOA_estimate)
    error_sum1=error_sum1+((DOA_estimate(i)-AL(i)))^2;
    end
    for i=length(DOA_estimate)+1:K
    error_sum1=error_sum1+((DOA_estimate(length(DOA_estimate))-AL(i)))^2;
    end
 end
 hold on
% figure
power=(real(Source_power_re)/max(real(Source_power_re)));
plot(Azimuth,real(power),'g-','Linewidth',0.5);
error1=error1+error_sum1/length(AL);
end
RMSE=sqrt(error1/test)  %root mean square error
hold on
scatter(AL(1),1,30,'filled','k');
scatter(AL(2),1,30,'filled','k');
scatter(AL(3),1,30,'filled','k');
box on
ylim([0 1])
xlim([-45 45])
xlabel('Azimuth','fontsize',12);
ylabel('unified spectrum','fontsize',12);
legend('CGDP-SBL')
function [source_power,Azimuth,source_id,iter]=CGDP_SBL(Y,Azimuth_grid,grid_interval,Aa,K_)
%% input:
%% Y: Array output data
%% Azimuth_grid: Space discontinue Azimuth angle grid
%% grid_interval: Interval of the adjacent angle grid 
%% Aa: Compeleted array steering matrix
%% K_: Source number,if K_ is not a prior,we let K_=M-1.

%% output:
%% Source_power:signal space power
%% Azimuth:All space angles, include estimated DOAs
%% iter:iter number
%% source_id: DOA index in Azimuth
%% Noise_variance: noise power
[M,L]=size(Y);
Ry=Y*Y'/L;
N=length(Azimuth_grid);
h=0.1; % parameter in (7)
%% Initialize
Noise_variance=0.1*(norm(Y))^2/(M*L); % noise variance
Beta=1/Noise_variance;
source_power_old = sum(abs(Aa'*Y),2)/(M*L); % signal power(gamma)
jmax=500;
xi=ones(N,1);
iter=0;
sigma_y=[];
gamma=[];
tol=0.001; % tolorance
while  (iter<jmax)
   iter=iter+1;  % number of iteraion     
   gamma=diag(source_power_old);
   B=Aa*gamma*Aa'; 
   sigma_y=1/Beta*eye(M)+B;    
   sigma_y_inv=inv(sigma_y);           
   Aa_sigma_y_inv=Aa'* sigma_y_inv;
   for i=1:N
      Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa(:,i);
   end
   sigma_x_ii=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  % eq.(22)
   Mu_x=gamma*Aa_sigma_y_inv*Y;  % eq.(11)
   Mu_x_norm= sum(abs(Mu_x).^2, 2);       
%% CGDP-SBL(eq.(15)):
   source_power=( (-2*L)+sqrt( (2*L)^2 +4*(Mu_x_norm+ L*real(sigma_x_ii)).*xi.^2) )./(xi.^2); 
   xi=( -h+sqrt( h^2+2*source_power*(h+2) ) )./(source_power);  % eq.(17)
   Beta=(M*L)/(((norm(Y-Aa*Mu_x,'fro'))^2)+L*trace(B-B*sigma_y_inv*B));  % eq.(18)
   if (norm((source_power-source_power_old),2)/norm(source_power_old,2)<tol) 
       break
   end
   source_power_old=source_power;
end
%%  refined DOA searching:
[~, peakindex1] = findpeaks(real(source_power),'sortstr','descend');
KP=min(length(peakindex1),K_);
Azimuth=Azimuth_grid;
source_id=sort( peakindex1(1:KP),'ascend');  
theta_r1=[];
Azimuth_refine1=[];
theta_refine1=[]; 
gamma_k=[];
sigma_yk_inv=[];
r_step=0.02;  % searching step
start2=tic;
Aa_re_=Aa;
gamma_re1=gamma;
sigma_y_re1=sigma_y;

for i=1:length(source_id)
    Aa_rek=[];
       ar=0;           
   if (source_power(source_id(i)-1)) < (source_power(source_id(i)+1))            
      Aa_rek= [Aa_re_(:, source_id(i):source_id(i)+1)];         
      gamma_k=diag(gamma_re1)';
      gamma_k= [gamma_k(source_id(i):source_id(i)+1)];
      gamma_k=diag(gamma_k);
      sigma_y_rek1=sigma_y_re1-Aa_rek* gamma_k*Aa_rek';    %
      sigma_yk_inv=inv( sigma_y_rek1); 
      qs=sigma_yk_inv*Ry*sigma_yk_inv*L;
           for agr= Azimuth(source_id(i)) :r_step:Azimuth(source_id(i))+grid_interval   % 
                 ar=ar+1;
                aa=[];ad=[];
                aa=exp(1j*pi*[0:M-1]'*sin(agr*pi/180));  % linear array steering vector  
                qk=aa'*qs*aa;zk=aa'*sigma_yk_inv*aa;
                source_power_rere(i,ar)= ( -xi(source_id(i))^2/2-L*zk+sqrt( ((xi(source_id(i))^2/2+L*zk))^2-xi(source_id(i))^2*(xi(source_id(i))^2/4+L*zk-qk)) )/(xi(source_id(i))^2/2*zk);    
                theta_r1(i,ar)= real( -L*log( ( 1+( source_power_rere(i,ar))*zk))+( qk)/( ( source_power_rere(i,ar))^(-1)+ zk )...
               -xi(source_id(i))^2/4*(source_power_rere(i,ar)));                                                         
           end
          
        [agr_max(i)]=find((theta_r1(i,:))==max(theta_r1(i,:)));   
         Azimuth_refine(i,:)= Azimuth(source_id(i)):r_step:Azimuth(source_id(i))+grid_interval;
        theta_refine(i)=Azimuth_refine(i,agr_max(i));  % eq.£¨21£© 
   else
   
      Aa_rek= [Aa_re_(:, source_id(i)-1:source_id(i))];        
      gamma_k=diag(gamma_re1)';
      gamma_k= [gamma_k(source_id(i)-1:source_id(i))];
      gamma_k=diag(gamma_k);
      sigma_y_rek1=sigma_y_re1-Aa_rek* gamma_k*Aa_rek';   
      sigma_yk_inv=inv( sigma_y_rek1);
      qs=sigma_yk_inv*Ry*sigma_yk_inv*L;
               aa=[];ad=[];
             for agr= Azimuth(source_id(i))-grid_interval:r_step:Azimuth(source_id(i))   
              ar=ar+1; 
              aa=exp(1j*pi*[0:M-1]'*sin(agr*pi/180)); 
              qk=aa'*qs*aa;zk=aa'*sigma_yk_inv*aa;
              source_power_rere(i,ar)= ( -xi(source_id(i))^2/2-L*zk+sqrt( ((xi(source_id(i))^2/2+L*zk))^2-xi(source_id(i))^2*(xi(source_id(i))^2/4+L*zk-qk)) )/(xi(source_id(i))^2/2*zk);  % eq£¨20£©  
              theta_r1(i,ar)= real( -L*log( ( 1+( source_power_rere(i,ar))*zk))+( qk)/( ( source_power_rere(i,ar))^(-1)+ zk )...
              -xi(source_id(i))^2/4*(source_power_rere(i,ar)));  
             end
      [agr_max(i)]=find((theta_r1(i,:))==max(theta_r1(i,:)));   
      Azimuth_refine(i,:)= Azimuth(source_id(i))-grid_interval:r_step:Azimuth(source_id(i));
      theta_refine(i)=Azimuth_refine(i,agr_max(i));  % eq.£¨21£© 
   end 
    Azimuth(source_id(i))=theta_refine(i);
end

end