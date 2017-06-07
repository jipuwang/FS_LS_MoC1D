% changes:
% maxIterate is changed to 500.
% N: number of angles is changed to 2
function [phi_j_new]=MoC1D_prob5_LS_sqrt(J,N)
%% Example also as optional param
if ~exist('J','var')
    J=5*2*2%*2*2%*2*2*2*2*2*2
    %J=5%*2%*2%*2%*2%*2%*2%*2%*2%*2%*2
end
% 1D MoC code
if ~exist('N','var')
    N=16;
end

Z=10;% 1D MoC code
%% input parameters
delta_z_j=ones(1,J)*Z/J;
% NOTE: CRIME HERE. THIS IS REQUIRING UNIFORM GRID ONLY!
delta_z=Z/J;
% The sum of the weight_n is 2. One CANNOT double the weight and cut the
% size of mu_n by half and trace only half of the angles.
% [mu_n,weight_n]=GSQuadSet(N);
% switched to new qaudrature set generator for higher order and faster
% performance.
[mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);

sigma_t=1.0;
sigma_s=0.5;

% set up boundary condition
% The below is isotropic incoming B.C.
psi_b1_n=ones(N,1)*1.0; % the first/negative half is not useful; n=N/2+1:N % mu>0
psi_b2_n=ones(N,1)*sqrt(Z+1); % the second/positive half is not useful; n=1:N/2 % mu<0
% for n=1:N/2
%     psi_b2_n(n)=1+Z*exp(mu_n(n));
% end
Q_n_j=ones(N,J)*1.0; % MMS source, angular density source
Q_n_j_hat=ones(N,J);
% overwrite with MMS source
% NOTE: ANOTHER CRIME
sigma_s_j=ones(1,J)*sigma_s;
sigma_t_j=ones(1,J)*sigma_t;
sigma_t_inv_j=1./sigma_t_j;

for j=1:J
    for n=1:N
        Q_n_j(n,j)=mu_n(n)*(sqrt(j*delta_z+1)-sqrt((j-1)*delta_z+1))/delta_z ...
            +(sigma_t-sigma_s)*(2/3*((j*delta_z+1)^1.5-((j-1)*delta_z+1)^1.5))/delta_z;

%             +(sigma_t-sigma_s)*(2/15*((j*delta_z+1)^1.5*(3*j*delta_z-2)-sqrt((j-1)*delta_z+1)*(3*(j-1)*delta_z*(j-1)*delta_z+(j-1)*delta_z-2)))/delta_z ...        
        Q_n_j_hat(n,j)=mu_n(n)*(1/3*((j*delta_z-2)*sqrt(j*delta_z+1)-((j-1)*delta_z-2)*sqrt((j-1)*delta_z+1)))/delta_z ...
            +(sigma_t-sigma_s)*(2/15*((j*delta_z+1)^1.5*(3*j*delta_z-2)-((j-1)*delta_z+1)^1.5*(3*(j-1)*delta_z-2)))/delta_z ...        
               -((j-1)+j)*delta_z*0.5*Q_n_j(n,j);
        % Simplification is possible to make the above delta_z*delta_z/12. 
    end % n
end % j

maxIterate=2000;
epsilon_phi=1e-14;

%% MoC1D Solver
% N rays to trace, each angle has only 1 ray, no ray-spacing
% n for each angle, and j for FSR region index
segLen_n_j=zeros(N,J);

for n=1:N
    for j=1:J
        segLen_n_j(n,j)=delta_z_j(j)/abs(mu_n(n));
    end
end

phi_j_old=ones(1,J)*1.0; % so the 1st dimension is consistently the angle. 
phi_j_old_hat=ones(1,J)*1.0; % so the 1st dimension is consistently the angle. 

Q_n_j_x=zeros(N,J); % these are actually angular quantities, already have 0.5's in them.
q_n_j=zeros(N,J);
q_n_j_sm=zeros(N,J);
% new quantities for linear source
Q_n_j_x_hat=zeros(N,J);
q_n_j_hat=zeros(N,J);
q_n_j_sm_hat=zeros(N,J);

for iIterate=1:maxIterate
    for j=1:J
        for n=1:N
            Q_n_j_x(n,j)=0.5*sigma_s_j(j)*phi_j_old(j)+Q_n_j(n,j);
            q_n_j(n,j)=Q_n_j_x(n,j);
            q_n_j_sm(n,j)=q_n_j(n,j);

            Q_n_j_x_hat(n,j)=0.5*sigma_s_j(j)*phi_j_old_hat(j)+Q_n_j_hat(n,j); 
            q_n_j_hat(n,j)=Q_n_j_x_hat(n,j)/(delta_z*delta_z/12);
            q_n_j_sm_hat(n,j)=q_n_j_hat(n,j)*(mu_n(n));   % NO ABS IS NEEDED!
        end
    end
%     phi_j_old_hat
%     q_n_j
%     q_n_j_hat

    phi_j_new=zeros(1,J);
    phi_j_new_hat=zeros(1,J);
    % ray tracing
    for n=N/2+1:N
        psi_in=psi_b1_n(n);
        for j=1:J
            tau_temp=sigma_t_j(j)*segLen_n_j(n,j);
            F1=1-exp(-tau_temp);
            F2=2*(tau_temp-F1)-tau_temp*F1;
            psi_out=psi_in+(q_n_j_sm(n,j)*sigma_t_inv_j(j)-psi_in)*F1...
                +(q_n_j_sm_hat(n,j)*0.5*sigma_t_inv_j(j)*sigma_t_inv_j(j))*F2;
            psi_avg=q_n_j_sm(n,j)*sigma_t_inv_j(j)+(psi_in-psi_out)/tau_temp;
            phi_j_new(j)=phi_j_new(j)+weight_n(n)*psi_avg;
            G1=1+tau_temp*0.5-(1+1/tau_temp)*F1;
            G2=2/3*tau_temp-(1+2/tau_temp)*G1;
            psi_hat=psi_in*0.5*segLen_n_j(n,j) ...
                + (q_n_j_sm(n,j)*sigma_t_inv_j(j)-psi_in)*G1*sigma_t_inv_j(j) ...
                + (q_n_j_sm_hat(n,j)*0.5*sigma_t_inv_j(j)*sigma_t_inv_j(j))*segLen_n_j(n,j)*G2; %changed here

            phi_j_new_hat(j)=phi_j_new_hat(j)+...
                weight_n(n)*(-delta_z_j(j)*0.5*psi_avg+abs(mu_n(n))*psi_hat);
            psi_in=psi_out;
        end
    end
    
    for n=1:N/2
        psi_in=psi_b2_n(n);
        for j=J:-1:1
            tau_temp=sigma_t_j(j)*segLen_n_j(n,j);
            F1=1-exp(-tau_temp);
            F2=2*(tau_temp-F1)-tau_temp*F1;
            psi_out=psi_in+(q_n_j_sm(n,j)*sigma_t_inv_j(j)-psi_in)*F1...
                +(q_n_j_sm_hat(n,j)*0.5*sigma_t_inv_j(j)*sigma_t_inv_j(j))*F2;
            psi_avg=q_n_j_sm(n,j)*sigma_t_inv_j(j)+(psi_in-psi_out)/tau_temp;
            phi_j_new(j)=phi_j_new(j)+weight_n(n)*psi_avg;
            G1=1+tau_temp*0.5-(1+1/tau_temp)*F1;
            G2=2/3*tau_temp-(1+2/tau_temp)*G1;
            psi_hat=psi_in*0.5*segLen_n_j(n,j) ...
                + (q_n_j_sm(n,j)*sigma_t_inv_j(j)-psi_in)*G1*sigma_t_inv_j(j) ...
                + (q_n_j_sm_hat(n,j)*0.5*sigma_t_inv_j(j)*sigma_t_inv_j(j))*segLen_n_j(n,j)*G2; %changed here

            phi_j_new_hat(j)=phi_j_new_hat(j)+...
                weight_n(n)*(+delta_z_j(j)*0.5*psi_avg-abs(mu_n(n))*psi_hat);
            psi_in=psi_out;
        end
    end
    
    % test for convergence
    error=norm(phi_j_new-phi_j_old);
    if error<epsilon_phi
        break;
    end
    phi_j_old=phi_j_new;
    phi_j_old_hat=phi_j_new_hat;
end
error
% phi_j_old=phi_j_new;
phi_j_new=phi_j_new';
display(iIterate);
% figure(19);
plot(phi_j_new,'*-')
% openvar('phi_j_new')
end
