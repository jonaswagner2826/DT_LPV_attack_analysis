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
openModel = false;
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
n = 4;%2; % States
m = 4; % Subsystems
p = 2;%1; % Inputs
q = 3;%1; % Outputs

% Polyotopic Matrix definitions
A(:,:,1) = [-0.80, 0.25; 0.25,-0.30]; B(:,:,1) = [ 1.90; 0.00];
A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00]; B(:,:,2) = [-1.00; 1.50];
A(:,:,3) = [-0.30, 0.65; 0.55, 0.10]; B(:,:,3) = [ 0.30;-2.00];
A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30]; B(:,:,4) = [-0.60; 0.00];

% A0 = [0,0;0,0]; B0 = [0;0];
% A1 = [-0.80, 0.25; 0.25,-0.30]; B1 = [ 1.90; 0.00];
% A2 = [ 0.30, 0.70; 0.70, 0.00]; B2 = [-1.00; 1.50];
% A3 = [-0.30, 0.65; 0.55, 0.10]; B3 = [ 0.30;-2.00];
% A4 = [ 0.55,-0.20;-0.40,-0.30]; B4 = [-0.60; 0.00];

a = -0.9;
b = 0.9;
A = a + (b-a) * rand(n,n,m);
B = a + (b-a) * rand(n,p,m);







% Output Eq
C = eye(q,n); D = zeros(q,p);
% C = [1, 0]; D = 0;

% Input
u0 = 0.5*rand(p,1);%0.5;
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

% DE Simulation Setup
Alpha_sym = sym('alpha',[m,1], 'real');
A_bar = zeros(n);
B_bar = zeros(n,p);
for i = 1:m
    A_bar = A_bar + Alpha_sym(i) * A(:,:,i);
    B_bar = B_bar + Alpha_sym(i) * B(:,:,i);
end

a_coeff_temp = charpoly(A_bar);
if size(a_coeff_temp,2) <= n
    a_coeff_temp = [zeros(1,n+1-size(a_coeff_temp,2)),a_coeff_temp];
end

temp = sym('temp');
for i = 1:q %iterating through outputs
    for j = 1:p %iterating through inputs
        clear b_coeff_temp
        b_coeff_temp = coeffs(C(i,:)*adjoint(temp*eye(n) - A_bar)*B_bar(:,j),temp);
        if size(b_coeff_temp,2) <= n
            b_coeff_temp = [zeros(1,n+1-size(b_coeff_temp,2)), b_coeff_temp];
        end
        a_coeff(:,i,j) = a_coeff_temp;
        b_coeff(:,i,j) = b_coeff_temp;
    end
end

nuVector = matlabFunction([a_coeff(1:n,:,:); b_coeff(1:n,:,:)], 'vars', {Alpha_sym});




% Attack (essentially just noise for this small system)
v0 = 0;
V = v0 * randn(q,k_max);



% Initial Conditions
x_0 = rand(n,1);%[0.25; -6.4];
x_hat_0 = rand(n,1);%[-1.25; 3.4];
alpha_hat_0 = [0.25; 0.25; 0.25; 0.25];



%% Rough Simulation of the System
% Sim Data
X = zeros(n,k_max);
Y = zeros(q, k_max);
Alpha = zeros(m,k_max);

% DE Data
X_hat_DE = zeros(n,k_max);
Alpha_hat_DE = zeros(m,k_max);
L_data_DE = zeros(n,q,m,k_max);
P_data_DE = zeros(2*n,2*n,k_max);
Phi = zeros(2*n,q,p,k_max);
Nu_hat = zeros(2*n,q,p,k_max);

% EKF Data
X_hat_EKF = zeros(n,k_max);
Alpha_hat_EKF = zeros(m,k_max);
P_data_EKF = zeros(n+m,n+m,k_max);


% Simulation Initialization
x = x_0;

% DE Initialization
x_hat_DE = x_hat_0;
alpha_hat_DE = alpha_hat_0;
L_DE = zeros(n,q,m);
% P_DE = randn(2*n); P_DE = P_DE*P_DE.'; %random ??? or this rand one:
P_DE = [1.7552    1.7244    0.7521    0.4716;
        1.7244    4.3679   -0.4500   -0.7000;
        0.7521   -0.4500    3.5146    2.5992;
        0.4716   -0.7000    2.5992    1.9602];
nu_hat = randn(2*n,q,p); %random ??? or this rand one:
% nu_hat = [0.3; -0.2; 0.1; -0.4];
% phi = randn(2*n,q,p); %random ??? or this rand one: %[1.2424; -1.0667; 0.9337; 0.3503];
phi = zeros(2*n,q,p);

gamma = 0.9;

% EKF Initialization
x_hat_EKF = x_hat_0;
alpha_hat_EKF = alpha_hat_0;
P_EKF = 1e5 * eye(n+m);% random large P...
Q_EKF = diag([zeros(1,n), 1*ones(1,m)]);
R_EKF = 0.01;

for k = K(1:100)
    %% Plant Simulation
    alpha = alpha_traj(k);
    u = U(:,k);
    x_old = x;
    x = 0;
    for i = 1:m
        x = x + alpha(i) * (A(:,:,i) * x_old + B(:,:,i) * u);
    end
    y = C*x + D*u + V(:,k);
    
    
    % Sim Data
    X(:,k) = x;
    Y(:,k) = y;
    Alpha(:,k) = alpha;

    %% State Estimation
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
    
    phi_old = phi;
    phi(2:n,:,:) = phi_old(1:n-1,:,:);
    phi(1,:,:) = y * ones(1,p);
    phi(n+2:2*n,:,:) = phi_old(n+1:2*n-1,:,:);
    phi(n+1,:,:) = (u * ones(1,q))';
    Phi(:,:,:,k) = phi;
    
    % CVX Feasability Problem on LMI
    disp('DE State Estimation Started')
    tol = 1e-6;
    cvx_clear
    cvx_begin sdp quiet
        variable P_cvx(n,n,m) symmetric
        variable F(n,q,m)
        variable G(n,n,m)
        variable zeta
        subject to
            for i = 1:m
                for j = 1:m
                    [G(:,:,i) + G(:,:,i)' - P_cvx(:,:,j), zeros(n),...
                        G(:,:,i)*A(:,:,i) + F(:,:,i)*C, G(:,:,i);
                     zeros(n), eye(n), eye(n), zeros(n);
                     A(:,:,i)'*G(:,:,i)' + C'*F(:,:,i)', eye(n),...
                        P_cvx(:,:,i), zeros(n);
                     G(:,:,i)', zeros(n), zeros(n), zeta*eye(n)
                     ]>= tol*eye(4*n);
                end
            end
            zeta >= 1;
    cvx_end

    % Gain Result Calc
    for i = 1:m
        L_DE(:,:,i) = G(:,:,i)\F(:,:,i);
    end
    
    % State Estimate Update
    x_hat_DE_old = x_hat_DE; % x_hat_{k-1}
    x_hat_DE = 0;
    for i = 1:m
        x_hat_DE = x_hat_DE + alpha_hat_DE(i)...
                        *(A(:,:,i)*x_hat_DE_old + B(:,:,i)*u)...
                    + L_DE(:,:,i)*(C*x_hat_DE_old - y);
    end
    L_data_DE(:,:,:,k) = L_DE;
    
    
    %% Parameter Estimation
               
    disp('DE Parameter Estimation Started')
    % attempting yalmip
    yalmip('clear')
    alpha_yalmip = sdpvar(m,1);
    nu_alpha = sdpvar(2*n,q,p);
    param_opt_sum = sdpvar(q,p,k);
    
    Constraints = [sum(alpha_yalmip) == 1];
    for i = 1:m
        Constraints  = [Constraints, alpha_yalmip(i) >= 0];
    end
    Constraints = [Constraints, nu_alpha == nuVector(alpha_yalmip)];
    
    for i = 1:q %iterating through outputs
        for j = 1:p %iterating through inputs
%             Constraints = [Constraints,...
%                            param_opt_sum(i,j,1) == norm(Y(i,1) - Phi(:,i,j,1) * nu_alpha(:,i,j))^2];    
            for t = 1:k
                Constraints = [Constraints, param_opt_sum(i,j,t) ...
                                    == gamma^(k-t) * norm(Y(i,t) ...
                                    - Phi(:,i,j,t)' * nu_alpha(:,i,j))^2];
            end
        end
    end
    
    Objective = sum(param_opt_sum,'all');
    
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
%     P_data_DE(:,:,k) = P_DE;
    Nu_hat(:,:,:,k) = nu_hat;
    
    
    %% EKF Method
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

%% Ploting
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
title('System Parameter Error')
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
title('System Parameter Error')
ax6 = subplot(3,2,6);
plot((C*X_hat_EKF)'-Y')
title('Output Error')

