% simulink Generation Scipt
% Jonas Wagner
% 2021-04-08
% 
% Important References:
% ----------------------
% https://www.mathworks.com/help/simulink/slref/add_block.html
% https://www.mathworks.com/help/simulink/programmatic-modeling.html
% https://www.mathworks.com/help/simulink/referencelist.html
% https://www.mathworks.com/help/simulink/slref/common-block-parameters.html
% https://www.mathworks.com/help/simulink/slref/block-specific-parameters.html

clear
close all

%% Settings
generateModel = false;
openModel = true;
simulateModel = false;
plotResults = false;

% Name of the simulink model
[cfolder,~,~] = fileparts(mfilename('fullpath'));
subfolder = ''; %include / at end of subfolder
fname = 'DT_LPV_sim';


%% System Definitions

%--------------------------------------------------------------------
% DT_LPV Model Definition
% (https://www.sciencedirect.com/science/article/pii/S2405896317313459)

% DT Settings
k_max = 100;
K = 1:k_max; % Simulation Time

% System Size
n = 2; % States
m = 4; % Subsystems
p = 1; % Inputs
q = 1; % Outputs

% Polyotopic Matrix definitions
A(:,:,1) = [-0.80, 0.25; 0.25,-0.30]; B(:,:,1) = [ 1.90; 0.00];
A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00]; B(:,:,2) = [-1.00; 1.50];
A(:,:,3) = [-0.30, 0.65; 0.55, 0.10]; B(:,:,3) = [ 0.30;-2.00];
A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30]; B(:,:,4) = [-0.60; 0.00];

A0 = [0,0;0,0]; B0 = [0;0];
A1 = [-0.80, 0.25; 0.25,-0.30]; B1 = [ 1.90; 0.00];
A2 = [ 0.30, 0.70; 0.70, 0.00]; B2 = [-1.00; 1.50];
A3 = [-0.30, 0.65; 0.55, 0.10]; B3 = [ 0.30;-2.00];
A4 = [ 0.55,-0.20;-0.40,-0.30]; B4 = [-0.60; 0.00];

% Output Eq
C = [1, 0]; D = 0;

% Input
u0 = 0.5;
T = 10;
U = zeros(1,k_max);
i = 1;
duty = 100;
for k = 1:k_max
    if i == T
        duty = 100*rand();
        i = 1;
    end
    U(:,k) = square((2*pi*k)/T, duty);
    i = i+1;
end
U = u0 * (U+1)/2;

% Attack (essentially just noise for this small system)
v0 = 0.5;
V = v0 * randn(q,k_max);



% Initial Conditions
x_0 = [0.25; -6.4];
x_hat_0 = [-1.25; 3.4];
alpha_hat_0 = [0.25; 0.25; 0.25; 0.25];



%% Rough Simulation of the System
% Sim Data
X = zeros(n,k_max);
Y = zeros(q, k_max);
Alpha = zeros(m,k_max);

% DE Data
X_hat_DE = zeros(n,k_max);
Alpha_hat_DE = zeros(m,k_max);
P_data_DE = zeros(2*n,2*n,k_max);
Phi = zeros(1,2*n,k_max);
Nu_hat = zeros(2*n,1,k_max);

% EKF Data
X_hat_EKF = zeros(n,k_max);
Alpha_hat_EKF = zeros(m,k_max);
P_data_EKF = zeros(n+m,n+m,k_max);


% Simulation Initialization
x = x_0;

% DE Initialization
x_hat_DE = x_hat_0;
alpha_hat_DE = alpha_hat_0;
P_DE = randn(2*n); P_DE = P_DE*P_DE.'; %random ??? or this rand one:
P_DE = [1.7552    1.7244    0.7521    0.4716;
        1.7244    4.3679   -0.4500   -0.7000;
        0.7521   -0.4500    3.5146    2.5992;
        0.4716   -0.7000    2.5992    1.9602];
% nu_hat = randn(2*n,1); %random ??? or this rand one:
nu_hat = [0.3; -0.2; 0.1; -0.4];
% phi_hat = randn(2*n,1); %random ??? or this rand one: %[1.2424; -1.0667; 0.9337; 0.3503];
phi = zeros(1,m);

gamma = 0.9;

% EKF Initialization
x_hat_EKF = x_hat_0;
alpha_hat_EKF = alpha_hat_0;
P_EKF = 1e5 * eye(n+m);% random large P...
Q_EKF = diag([zeros(1,n), 100*ones(1,m)]);
R_EKF = 0.01;

for k = K(1:100)
    % Plant Simulation
    alpha = alpha_traj(k);
    u = U(k);
    x_old = x;
    x = 0;
    for i = 1:m
        x = x + alpha(i) * (A(:,:,i) * x_old + B(:,i) * u);
    end
    y = C*x + D*u + V(:,k);
    
    
    % Sim Data
    X(:,k) = x;
    Y(:,k) = y;
    Alpha(:,k) = alpha;

    % DE Method
    if k >= 3
    disp(['DE Method: k = ', num2str(k)])
    P_DE_old = P_DE;
    P_DE = (1/gamma)*P_data_DE(:,:,k-1) - (1/gamma)*(phi'...
        * (gamma + phi * P_data_DE(:,:,k-1) * phi')^(-1)*(phi * P_data_DE(:,:,k-1)));
    if any(eig(inv(P_DE)) < 0)
        warning('INV(P_DE) not PSD')
        
    end
    
    end
%     phi = [fliplr(Y(:,(k-2):(k-1))), fliplr(U((k-2):(k-1)))];
    phi_old = phi;
    phi(2:n) = phi_old(1:n-1);
    phi(1) = y;
    phi(n+2:2*n) = phi_old(n+1:2*n-1);
    phi(n+1) = u;
    Phi(:,:,k) = phi;
    
    
    [x_hat_DE] = est_DE(x_hat_DE, alpha_hat_DE, y, u,...
                                           P_DE, phi, nu_hat,...
                                           gamma, A, B, C);
    
                      
    % Param. Opt. Problem
    disp('DE Parameter Estimation Started')
    % attempting yalmip
    yalmip('clear')
    alpha_yalmip = sdpvar(m,1);
    nu_alpha = sdpvar(2*n,1);
    param_opt_sum = sdpvar(k,1);
    
    Constraints = [sum(alpha_yalmip) == 1];
    Constraints = [Constraints, nu_alpha == nuVector(alpha_yalmip)];
    
    Constraints = [Constraints, nu_alpha == [
    1;
    -((3*alpha_yalmip(2))/10 - (11*alpha_yalmip(1))/10 - alpha_yalmip(3)/5 + alpha_yalmip(4)/4 + 1);
    0;
    (alpha_yalmip(2) - (19*alpha_yalmip(1))/10 - (3*alpha_yalmip(3))/10 + (3*alpha_yalmip(4))/5 + 1)
    ]];
    
    
    for i = 1:m
        Constraints  = [Constraints, alpha_yalmip(i) >= 0];
    end
    Constraints = [Constraints, param_opt_sum(1) == norm(Y(:,1) - Phi(1) * nu_alpha)^2];
        
    for i = 2:k
        Constraints = [Constraints, param_opt_sum(i) ...
            == gamma * param_opt_sum(k-1) + norm(Y(:,k) - Phi(k) * nu_alpha)^2];
    end
    
    Objective = sum(param_opt_sum);
    
    options = sdpsettings();
    
    sol = optimize(Constraints,Objective,options);

    if sol.problem == 0
         % Extract and display value
         alpha_hat_DE = value(alpha_yalmip);
         nu_hat = value(nu_alpha);
    else
         disp('Hmm, something went wrong!');
         sol.info
         yalmiperror(sol.problem)
         error('issue here.....')
    end
    
    % DE Data
    X_hat_DE(:,k) = x_hat_DE;
    Alpha_hat_DE(:,k) = alpha_hat_DE;
    P_data_DE(:,:,k) = P_DE;
    Nu_hat(:,:,k) = nu_hat;
    
    % EKF Method
    disp(['EKF Method: k = ', num2str(k)])
    [x_hat_EKF, alpha_hat_EKF, P_EKF] = est_EKF(x_hat_EKF,...
                                                alpha_hat_EKF,...
                                                P_EKF, y, u,...
                                                A, B, C, Q_EKF, R_EKF);
    % EKF Data
    X_hat_EKF(:,k) = x_hat_EKF;
    Alpha_hat_EKF(:,k) = alpha_hat_EKF;
    P_data_EKF(:,:,k) = P_EKF;
    
end


X = X(:,1:k);
Alpha = Alpha(:,1:k);
Y = Y(:,1:k);
X_hat_DE = X_hat_DE(:,1:k);
Alpha_hat_DE = Alpha_hat_DE(:,1:k);
X_hat_EKF = X_hat_EKF(:,1:k);
Alpha_hat_EKF = Alpha_hat_EKF(:,1:k);

figure('position', [0,0,500,750])
sgtitle(['Attack Power: ',num2str(v0)])
subplot(3,1,1)
plot(X')
title('System States')
subplot(3,1,2)
plot(Alpha')
title('System Parameters')
subplot(3,1,3)
plot(Y')
title('System Output')


figure('position', [0,0,1000,750])
sgtitle(['DE: Attack Power = ',num2str(v0)])
ax1 = subplot(3,2,1);
plot(X_hat_DE')
title('System State Estimates')
ax2 = subplot(3,2,3);
plot(Alpha_hat_DE')
title('System Parameter Estimates')
ax3 = subplot(3,2,5);
plot((C*X_hat_DE)')
title('Output Estimates')
ax4 = subplot(3,2,2);
plot(X_hat_DE'-X')
title('System State Error')
ax5 = subplot(3,2,4);
plot(Alpha_hat_DE'- Alpha')
title('System Parameter Estimates')
ax6 = subplot(3,2,6);
plot((C*X_hat_DE)'-Y')
title('Output Error')

figure('position', [0,0,1000,750])
sgtitle(['EKF: ','Attack Power = ',num2str(v0)])
ax1 = subplot(3,2,1);
plot(X_hat_EKF')
title('System State Estimates')
ax2 = subplot(3,2,3);
plot(Alpha_hat_EKF')
title('System Parameter Estimates')
ax3 = subplot(3,2,5);
plot((C*X_hat_EKF)')
title('Output Estimates')
ax4 = subplot(3,2,2);
plot(X_hat_EKF'-X')
title('System State Error')
ax5 = subplot(3,2,4);
plot(Alpha_hat_EKF'- Alpha')
title('System Parameter Estimates')
ax6 = subplot(3,2,6);
plot((C*X_hat_EKF)'-Y')
title('Output Error')

