%% CGDP-SBL2
clear all
close all
%
M=10;          % sensor number
T=40;          % snapshot number
SNR=-2;         % SNR
tt=1:T;
L=length(tt);
error=0;error1=0;
test=200; % test number
timetol_tol=0;
iter_tol=0;% iter number
for i=1:test  % independent simulation
    grid_interval=2;  % Azimuth grid_interval
    Azimuth_grid=-90:grid_interval:90;  % all hypothetical DOA angles
    %% 1.random DOAs£º
    u=unifrnd(-grid_interval/2,grid_interval/2);
    AL=[ -7+u 3+u 10+u ];
    %% 2.fixed DOAs
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
    
    %%
    tstart1 = tic;
    [Source_power_re,Azimuth,source_id_,j1]=CGDP_SBL2(Y,Azimuth_grid,grid_interval,Aa,K);
% [Source_power_re,Azimuth,source_id_,j1]=CGDP_SBL2(Y,Azimuth_grid,grid_interval,Aa,M-1);% if K is unknown
    time1=toc(tstart1);
    timetol_tol=timetol_tol+time1;  %total time
    [~, S_ID] = findpeaks((Source_power_re),'sortstr','descend');
    KP=min(K,length(S_ID));
    source_ID=S_ID(1:KP);
    source_ID=sort(source_ID,'ascend');
    DOA_estimate=Azimuth(source_ID(1:KP));
    
    iter_tol=iter_tol+j1; % iteration number
    error_sum1=0;
    if length( DOA_estimate)==K
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
    hold on
    scatter(AL(1),1,20,'filled','k');
    scatter(AL(2),1,20,'filled','k');
    scatter(AL(3),1,20,'filled','k');
end
iter_mean=iter_tol/test % mean iter number
RMSE=sqrt(error1/test)  %root mean square error

box on
ylim([0 1])
xlim([-45 45])
xlabel('Azimuth','fontsize',12);
ylabel('unified spectrum','fontsize',12);
legend('CGDP-SBL')
function [source_power,Azimuth,source_id,iter]=CGDP_SBL2(Y,Azimuth_grid,grid_interval,Aa,K)
%% input:
%% Y: Array output data
%% Azimuth_grid: Space discontinue Azimuth angle grid
%% grid_interval: Interval of the adjacent angle grid
%% Aa: Compeleted array steering matrix

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
Noise_variance=10^(-2)*(norm(Y))^2/(M*L); % noise variance
Beta=1/Noise_variance;
source_power_old = sum(abs(Aa'*Y),2)/(M*L); % signal power(gamma)
jmax=500;
xi=ones(N,1);
% xi=( -h+sqrt( h^2+2*source_power_old*(h+2) ) )./(source_power_old);  % also can be used 
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
    %% CGDP-SBL2:
    source_power=( 2*L*( -1+(real(sigma_x_ii))./source_power_old)+1+2*sqrt(( L*(-1+(real(sigma_x_ii))./source_power_old) +1/2).^(2)+xi.^2.*Mu_x_norm)) ./(xi.^2); % eq.(15£©
    xi=( -h+sqrt( h^2+2*source_power*(h+2) ) )./(source_power);  % eq.(16)
    source_id=[];
    source_power1=[];source_power2=[];source_power3=[];
    source_power1 = source_power(1:end-2);
    source_power2 = source_power(2:end-1);
    source_power3 = source_power(3:end);
    IX=[]; PeakLoc0=[];SortIdx=[];SigLoc=[];source_id=[];SigLoc_=[];source_id_=[];
    IX = find((source_power2 > source_power1) + (source_power2 >= source_power3) > 1);
    PeakLoc0 = IX+1;   %PeakLoc0 is index of spectrum
    [~,SortIdx] = sort(source_power(PeakLoc0),'descend');%descend the spectrum£¬SortIdx
    KP = min(length(PeakLoc0),K);
    SigLoc = PeakLoc0(SortIdx(1:KP));
    source_id_=sort(SigLoc,'ascend');
    signal_id=source_id_;
    KP_=length(signal_id);
    
    Ac=[Aa(:,signal_id)];
    P=[];
    P=Ac*(Ac'*Ac)^(-1)*Ac';
    Beta=1/abs(trace((eye(M)-P)*Ry)/(M-KP_));% eq.(18)
    
    
    if (norm((source_power-source_power_old),2)/norm(source_power_old,2)<tol)
        break
    end
    source_power_old=source_power;
end
%%  refined DOA searching:
[~, peakindex1] = findpeaks(real(source_power),'sortstr','descend');
KP=min(length(peakindex1),K);
Azimuth=Azimuth_grid;
source_id=sort( peakindex1(1:KP),'ascend');
theta_r1=[];
Azimuth_refine1=[];
theta_refine1=[];
gamma_k=[];
sigma_yk_inv=[];
r_step=0.05;  % searching step
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