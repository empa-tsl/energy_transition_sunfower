tic
clear variables
close all
clc

%% simulation setup

t_L = 30; %a (lifetime)
EROI = 20; %- (energy return on investment)

dt = .01; %a (time step)
time_delay = 4; %a from 1.1.2018; assumed that fossil engine runs constant in this time.

alpha = [0 0.5 1]; % - (continuous replacement of fossil in the economy)
n_alpha = size(alpha,2);

P_demand = 6; %TW (in 2018 without renewable fraction, which doesn't need to be replaced)
dot_m_CO2_current = 35; %Gt/a
%beta_step = 0.001;
beta = 0.4;
n_beta = size(beta,1);

delta_t_i = [0 0.01 0.02 0.04 0.1]; %a independence time for battery
n_delta_t_i = size(delta_t_i,2);

EI_storage = [60 80 400]; %- energy intensity: embodied electrical energy / storage capactity
eta_turnaround = [0.16 0.7 0.85];
n_storage_technology = size(EI_storage,2);


% required oversize as function of storage losses
phi_storage = [0.2 0.4 0.6]; %- fraction of harvested energy needing storage
%omega = [2]; % oversizing factor for the solar engine to meet energy losses in required storage and transmission.
n_phi = size(phi_storage,2);




%% simulation

P_invest = beta.*P_demand; %TW
t = 0:dt:100; %a
n_t = size(t,2);

P = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
EPBT = zeros(1,n_alpha,1,n_delta_t_i,n_storage_technology);
tau = zeros(1,n_alpha,1,n_delta_t_i,n_storage_technology);

for p=1:n_phi
    for i=1:n_delta_t_i
        for tech=1:n_storage_technology
            for k=1:n_alpha
            
            
                if i==1
                    EPBT(1,k,p,i,tech) = (t_L/EROI);
                else
                    EPBT(1,k,p,i,tech) = (t_L/EROI*(1+phi_storage(1,p)*alpha(1,k)*(1/eta_turnaround(1,tech)-1)) + delta_t_i(1,i) * EI_storage(1,tech)*phi_storage(1,p)/sqrt(eta_turnaround(1,tech))); % a (energy payback time)
                end
                
                tau(1,k,p,i,tech) = EPBT(1,k,p,i,tech)./((1-alpha(1,k))); % adjusted doubling time constant
                
            end
        end
    end
end

omega=zeros(n_phi,n_delta_t_i,n_storage_technology);

for i=1:n_delta_t_i
    for tech=1:n_storage_technology
        for p=1:n_phi
            if delta_t_i(1,i) == 0
                omega(p,i,tech) = 1;
            else
                omega(p,i,tech) = 1+phi_storage(1,p)*(1/eta_turnaround(1,tech)-1);
            end
                       
            P(1,:,p,i,tech) = 0;
            P(2,:,p,i,tech) = P_invest * t(1,2) / EPBT(1,1,p,i,tech);
            for j=3:n_t
                for k=1:n_alpha
                  
                    P(j,k,p,i,tech) = min([P(j-1,k,p,i,tech) + P(2,k,p,i,tech) .* 2.^(t(1,j-1)./tau(1,k,p,i,tech)),P_demand * omega(p,i,tech)]);
                end
            end
        end
    end
end

P_rep = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);

for i=1:n_delta_t_i
    for tech=1:n_storage_technology
        for p=1:n_phi
            for j=1:n_t
                for k=1:n_alpha
                    if alpha(1,k) * P (j,k,p,i,tech) >=  P_demand
                        P_rep(j,k,p,i,tech) = P_demand;
                    else
                        P_rep(j,k,p,i,tech) = alpha(1,k) .* P (j,k,p,i,tech);
                    end
                end
            end
        end
    end
end

% Calculating transition time and CO2 emissions
t_transition = t(1,end)*ones(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
dot_m_CO2 = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
m_CO2 = zeros(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
for i=1:n_delta_t_i
    for tech=1:n_storage_technology
        for p=1:n_phi
            for k=1:n_alpha
                for j=1:n_t
                    if P(j,k,p,i,tech) >= omega(p,i,tech) * P_demand
                        t_transition(1,k,p,i,tech) = t(1,j);
                        break
                    end
                
                    if t(1,j) > t_transition(1,k,p,i,tech)
                        dot_m_CO2(j,k,p,i,tech) = 0;
                    else
                        dot_m_CO2(j,k,p,i,tech) = dot_m_CO2_current * (1 + beta - P_rep(j,k,p,i,tech)/ P_demand);
                    end
                end             
                m_CO2(1,k,p,i,tech) = sum(dot_m_CO2(:,k,p,i,tech) .* dt,1) + dot_m_CO2_current/0.8246 * time_delay; % also land use and other emissions stay constant during time delay
                if t_transition(1,k,p,i,tech) == t(1,end)
                    m_CO2(1,k,p,i,tech) = 5000 + m_CO2(1,k,p,i,tech); % if transition cannot be completed, cumulative emissions go to infinity for transition
                end
            end
        end
    end
end


toc