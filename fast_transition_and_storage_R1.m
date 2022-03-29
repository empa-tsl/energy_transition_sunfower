tic
clear variables
close all
clc

%% simulation setup

t_L = 30; %a (lifetime)
EROI = 20; %- (energy return on investment)

EPBT_PV0 = t_L/EROI;
LR_PV = 0.128; %- learning rate of PV for embodied energy
P_0=0.07; %TW installed solar output capacity at the beginning of transition

dt = .01; %a (time step)
time_delay = 3; %a from 1.1.2020; (for SR1.5: 1.1.2018) assumed that fossil engine runs constant in this time.

alpha_step=0.5;
alpha = [0:alpha_step:1]; % - (continuous replacement of fossil in the economy)
n_alpha = size(alpha,2);

P_demand = 6; %TW (in 2018 without renewable fraction, which doesn't need to be replaced)
dot_m_CO2_current = 35; %Gt/a
beta = 0.4; % fossil investment
n_beta = size(beta,1);

delta_t_i_step = 0.001;
delta_t_i = [0 0.002:delta_t_i_step:0.1]; %a independence time for battery
n_delta_t_i = size(delta_t_i,2);

EI_storage0 = [63 85 259]; %- energy intensity: embodied electrical energy / storage capactity

eta_out = [0.5 0.9 0.97]; %- input efficiency
eta_in = [0.43 0.84 0.97]; % output efficiency
eta_roundtrip = eta_out .* eta_in; %- storage round-trip efficiency
n_storage_technology = size(EI_storage0,2);

LR_storage = [0.05 0 0.1]; %-learning rate of energy storage technologies for embodied energy

E_storage0 = [5.8*10^(-7) 0.000133 1.13*10^(-4)]; % TWa installed storage capacity at beginning of transition



phi_step=0.1;
phi_storage = [0.2:phi_step:0.8]; %- fraction of harvested energy needing storage
n_phi = size(phi_storage,2);


%%
t = 0:dt:100; %a
n_t = size(t,2);


%% simulation

P_invest = beta.*P_demand; %TW


P_L = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology); %TW power supply capacity, with learning
P_NL = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology); %TW power supply capacity, without learning
EPBT = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
tau = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
EPBT_PV = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
E_storage_demand = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
EI = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);

omega=zeros(n_phi,n_delta_t_i,n_storage_technology);

for tech=1:n_storage_technology
    for i=1:n_delta_t_i
        for p=1:n_phi
            if delta_t_i(1,i) == 0
                omega(p,i,tech) = 1;
            else
                omega(p,i,tech) = 1+phi_storage(1,p)*(1/eta_roundtrip(1,tech)-1);
            end
            for k=1:n_alpha
                % t=0
                EPBT_PV(1,k,p,i,tech) = EPBT_PV0;
                EPBT_PV(2,k,p,i,tech) = EPBT_PV0;
                EI(1,k,p,i,tech) = EI_storage0(1,tech);
                EI(2,k,p,i,tech) = EI_storage0(1,tech);
                if i==1
                    EPBT(1,k,p,i,tech) = EPBT_PV0;
                else
                    EPBT(1,k,p,i,tech) = (EPBT_PV0*(1+phi_storage(1,p)*alpha(1,k)*(1/eta_roundtrip(1,tech)-1)) + delta_t_i(1,i) * EI_storage0(1,tech)/eta_out(1,tech)); % a (energy payback time)
                end
              
                P_L(1,k,p,i,tech) = 0;
                P_NL(1,k,p,i,tech) = 0;
                % t=dt
                EPBT(2,:,p,i,tech)=EPBT(1,:,p,i,tech);
                P_L(2,k,p,i,tech) = P_invest * dt ./ EPBT(1,k,p,i,tech);
                P_NL(2,k,p,i,tech) = P_invest * dt ./ EPBT(1,k,p,i,tech);
                % t>dt
                for j=3:n_t
                    EPBT_PV(j,k,p,i,tech) = EPBT_PV0*((P_L(j-1,k,p,i,tech)+P_0)/P_0)^(log10(1-LR_PV)/log10(2));
                    E_storage_demand(j,k,p,i,tech) = P_L(j-1,k,p,i,tech) /(1+ phi_storage(1,p)*(1/eta_roundtrip(1,tech)-1))* delta_t_i(1,i)/(eta_out(1,tech));
                    EI(j,k,p,i,tech) = EI_storage0(1,tech)*((E_storage_demand(j,k,p,i,tech)+E_storage0(1,tech))/E_storage0(1,tech))^(log10(1-LR_storage(1,tech))/log10(2));
                    if i==1
                        EPBT(j,k,p,i,tech) = EPBT_PV(j,k,p,i,tech);
                    else
                        EPBT(j,k,p,i,tech) = (EPBT_PV(j,k,p,i,tech)*(1+phi_storage(1,p)*alpha(1,k)*(1/eta_roundtrip(1,tech)-1)) + delta_t_i(1,i) * EI(j,k,p,i,tech)/(eta_out(1,tech)));%*phi_storage(1,p))); % a (energy payback time)
                        
                    end
                    tau(j,k,p,i,tech) = EPBT(j,k,p,i,tech)./((1-alpha(1,k))); % adjusted doubling time constant

                    P_L(j,k,p,i,tech) = min([P_L(j-1,k,p,i,tech) + (P_invest + P_L(j-1,k,p,i,tech)*(1-alpha(1,k))) .* dt/EPBT(j,k,p,i,tech),P_demand * omega(p,i,tech)]);
                    P_NL(j,k,p,i,tech) = min([P_NL(j-1,k,p,i,tech) + (P_invest + P_NL(j-1,k,p,i,tech)*(1-alpha(1,k))) .* dt/EPBT(1,k,p,i,tech),P_demand * omega(p,i,tech)]);
                end
            end
        end
    end
end


%%
P_rep = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
P_rep_NL = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);

for i=1:n_delta_t_i
    for tech=1:n_storage_technology
        for p=1:n_phi
            for j=1:n_t
                for k=1:n_alpha
                    if alpha(1,k) * P_L (j,k,p,i,tech) >=  P_demand
                        P_rep(j,k,p,i,tech) = P_demand;
                    else
                        P_rep(j,k,p,i,tech) = alpha(1,k) .* P_L (j,k,p,i,tech);
                    end
                    if alpha(1,k) * P_NL (j,k,p,i,tech) >=  P_demand
                        P_rep_NL(j,k,p,i,tech) = P_demand;
                    else
                        P_rep_NL(j,k,p,i,tech) = alpha(1,k) .* P_NL (j,k,p,i,tech);
                    end
                end
            end
        end
    end
end

%% Calculating transition time and CO2 emissions
t_transition = t(1,end)*ones(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
dot_m_CO2 = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
m_CO2 = zeros(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);

t_transition_NL = t(1,end)*ones(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
dot_m_CO2_NL = zeros(n_t,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
m_CO2_NL = zeros(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
for i=1:n_delta_t_i
    for tech=1:n_storage_technology
        for p=1:n_phi
            for k=1:n_alpha
                for j=1:n_t
                    if P_L(j,k,p,i,tech) >= omega(p,i,tech) * P_demand
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
                    m_CO2(1,k,p,i,tech) = 1000 + m_CO2(1,k,p,i,tech); % if transition cannot be completed, cumulative emissions go to infinity for transition
                end
                % no learning
                for j=1:n_t
                    if P_NL(j,k,p,i,tech) >= omega(p,i,tech) * P_demand
                        t_transition_NL(1,k,p,i,tech) = t(1,j);
                        break
                    end
                
                    if t(1,j) > t_transition_NL(1,k,p,i,tech)
                        dot_m_CO2_NL(j,k,p,i,tech) = 0;
                    else
                        dot_m_CO2_NL(j,k,p,i,tech) = dot_m_CO2_current * (1 + beta - P_rep_NL(j,k,p,i,tech)/ P_demand);
                    end
                end             
                m_CO2_NL(1,k,p,i,tech) = sum(dot_m_CO2_NL(:,k,p,i,tech) .* dt,1) + dot_m_CO2_current/0.8246 * time_delay; % also land use and other emissions stay constant during time delay
                if t_transition_NL(1,k,p,i,tech) == t(1,end)
                    m_CO2_NL(1,k,p,i,tech) = 1000 + m_CO2_NL(1,k,p,i,tech); % if transition cannot be completed, cumulative emissions go to infinity for transition
                end
            end
        end
    end
end

%% CDF for 1.5°

load('CDF_1_5.mat')

%% Probability of violation

P_v_1_5 = zeros(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
P_v_1_5_NL = zeros(1,n_alpha,n_phi,n_delta_t_i,n_storage_technology);
for i=1:n_delta_t_i
    for tech=1:n_storage_technology
        for p=1:n_phi
            for k=1:n_alpha
                P_v_1_5(1,k,p,i,tech) = CDF_m_CO2_remaining_1_5_AR(1,round(m_CO2(1,k,p,i,tech)));
                P_v_1_5_NL(1,k,p,i,tech) = CDF_m_CO2_remaining_1_5_AR(1,round(m_CO2_NL(1,k,p,i,tech)));
            end
        end
    end
end


%%
toc
