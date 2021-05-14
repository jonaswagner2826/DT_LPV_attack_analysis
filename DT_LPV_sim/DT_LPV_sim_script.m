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
k_max = 1000;
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
T = 10;
U = zeros(size(K,1));
for i = floor(size(K,1)/T)
    for j = 1:T
        U(T*i + j) = square((2*pi*(T*i+j))/T);
    end
end

% Scedule Parameter
syms k_sym
g = piecewise(k_sym <  500, [0.50; 0.30; 0.20; 0.00],...
                     k_sym >= 500, [0.35; 0.40; 0.10; 0.15]);


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
Phi = zeros(2*n,k_max);
Nu_hat = zeros(2*n,k_max);

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
% P_DE = [1.7552    1.7244    0.7521    0.4716;
%         1.7244    4.3679   -0.4500   -0.7000;
%         0.7521   -0.4500    3.5146    2.5992;
%         0.4716   -0.7000    2.5992    1.9602];
% nu_hat = randn(2*n,1); %random ??? or this rand one:
nu_hat = [0.3; -0.2; 0.1; -0.4];
% phi_hat = randn(2*n,1); %random ??? or this rand one: %[1.2424; -1.0667; 0.9337; 0.3503];
phi = zeros(1,m);

gamma = 0.9;

% Alpha_sym = sym('alpha',[m,1]);
% A_bar = zeros(n);
% B_bar = zeros(n,p);
% for i = 1:m
%     A_bar = A_bar + Alpha_sym(i) * A(:,:,i);
%     B_bar = B_bar + Alpha_sym(i) * B(:,:,i);
% end
% a_coeff = charpoly(A_bar);
% temp = sym('temp');
% b_coeff = charpoly(C*adjoint(temp*eye(n) - A_bar)*B_bar);
% if size(a_coeff,2) <= n
%     a_coeff = [zeros(1,n+1-size(a_coeff,2)),a_coeff];
% end
% if size(b_coeff,2) <= n
%     b_coeff = [zeros(1,n+1-size(b_coeff,2)),b_coeff];
% end
% nu_alpha = matlabFunction([fliplr(a_coeff(1:n)), fliplr(b_coeff(1:n))]');


% EKF Initialization
x_hat_EKF = x_hat_0;
alpha_hat_EKF = alpha_hat_0;
P_EKF = 1e5 * eye(n+m);% random large P...
Q_EKF = diag([zeros(1,n), 100*ones(1,m)]);
R_EKF = 0.01;

% IMM Initialization


for k = K(1:100)
    % Plant Simulation
    alpha = alpha_traj(k);
    u = U(k);
    x_old = x;
    x = 0;
    for i = 1:m
        x = x + alpha(i) * (A(:,:,i) * x_old + B(:,i) * u);
    end
    y = C*x + D*u;
    
    % Sim Data
    X(:,k) = x;
    Y(:,k) = y;
    Alpha(:,k) = alpha;

    % DE Method
    if k >= 3
    disp(['DE Method: k = ', num2str(k)])
    phi = [fliplr(Y(:,(k-2):(k-1))), fliplr(U((k-2):(k-1)))];
    P_DE = (1/gamma)*P_data_DE(:,:,k-1) - (1/gamma)*(phi'...
        * (gamma + phi * P_data_DE(:,:,k-1) * phi')^(-1)*(phi * P_data_DE(:,:,k-1)));
    [x_hat_DE, alpha_hat_DE, nu_hat]...
                                = est_DE(x_hat_DE, alpha_hat_DE, y, u,...
                                           P_DE, phi, nu_hat,...
                                           gamma, A, B, C);
    end
    % DE Data
    X_hat_DE(:,k) = x_hat_DE;
    Alpha_hat_DE(:,k) = alpha_hat_DE;
    P_data_DE(:,:,k) = P_DE;
    Phi(:,k) = phi';
    Nu_hat(:,k) = nu_hat;
    
    % EKF Method
    disp(['EKF Method: k = ', num2str(k)])
    [x_hat_EKF, alpha_hat_EKF, P_EKF] = est_EKF(x_hat_EKF,...
                                                alpha_hat_EKF,...
                                                P_EKF, y, u,...
                                                A, B, C, Q_EKF, R_EKF);
    % IMM Method
    
                                            
                                            
    
    % Save Data to Arrays    

    
    
    % EKF Data
    X_hat_EKF(:,k) = x_hat_EKF;
    Alpha_hat_EKF(:,k) = alpha_hat_EKF;
    P_data_EKF(:,:,k) = P_EKF;
    
end

figure()
sgtitle('DE')
subplot(2,1,1)
plot(X_hat_DE'-X(:,1:size(X_hat_DE,2))')
subplot(2,1,2)
plot(Alpha(:,1:size(Alpha_hat_DE,2))'-Alpha_hat_DE')%


figure()
sgtitle('EKF')
subplot(2,1,1)
plot(X_hat_EKF'-X')
subplot(2,1,2)
plot(Alpha_hat_EKF'-Alpha')

% if generateModel
% %% Simulink Creation In Code
% % Simulink Settings ----------------------
% % Get the current configuration
% cfg = Simulink.fileGenControl('getConfig');
% % Changes Code Save Location
% cfg.CacheFolder = [pwd, '\', subfolder];
% cfg.CodeGenFolder = [pwd, '\', subfolder];
% cfg.CodeGenFolderStructure = 'TargetEnvironmentSubfolder';
% % Apply new Config
% Simulink.fileGenControl('setConfig', 'config', cfg, 'keepPreviousPath',true, 'createDir',true);
% 
% % Check if the file already exists and delete it if it does
% if exist(fname,'file') == 4
%     % If it does then check whether it's open
%     if bdIsLoaded(fname)
%         % If it is then close it (without saving!)
%         close_system(fname,0)
%     end
%     % delete the file
%     delete([fname,'.slx']);
% end
% 
% % Create Simulink Model
% new_system;
% 
% 
% 
% 
% 
% 
% % 
% % % Create Simiple Input
% % add_block('simulink/Sources/In1', [gcs, '/In']);
% % 
% % % Create Sum block
% % add_block('simulink/Commonly Used Blocks/Sum', [gcs, '/Sum'],...
% %     'inputs', '|+-');
% % 
% % % Connections (Simple Input output 1 to Sum Block input 1)
% % add_line(gcs, 'In/1', 'Sum/1'); 
% % 
% % % Saturation Block
% % add_block('simulink/Commonly Used Blocks/Saturation', [gcs, '/Saturation'], ...
% %     'LowerLimit','a', ... % block parameters
% %     'UpperLimit','b');
% % add_line(gcs, 'Sum/1', 'Saturation/1');
% % 
% % % State-Space System
% % add_block('cstblocks/LTI System', [gcs, '/LTI_sys'],...
% %     'sys','lti_sys',...
% %     'IC', 'x0');
% % add_line(gcs, 'Saturation/1', 'LTI_sys/1');
% % 
% % % Controller (just a feedback gain)
% % add_block('simulink/Commonly Used Blocks/Gain', [gcs, '/Controller'],...
% %     'Gain', 'K',...
% %     'Multiplication', 'Matrix(K*u)',...
% %     'Orientation', 'left');
% % add_line(gcs, 'LTI_sys/1', 'Controller/1');
% % add_line(gcs, 'Controller/1', 'Sum/2');
% % 
% % % Create Simple Scope/Output
% % add_block('simulink/Sinks/Out1', [gcs, '/Out']);
% % add_line(gcs, 'LTI_sys/1', 'Out/1');
% 
% 
% 
% %--------------------------------------------------------------------
% % LPV System
% 
% % DT_LPV_sys
% add_block('cstblocks/Linear Parameter Varying/DT LPV System',...
%     [gcs, '/DT_LPV_sys'])
% get_param(gcs, 'DT_LPV_sys/State-space array')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Auto Arrange
% Simulink.BlockDiagram.arrangeSystem(gcs) %Auto Arrange
% 
% %% Save and Open System
% save_system(gcs,[cfolder, '\', subfolder, fname]);
% print(['-s', gcs], '-dpng',... % Print model to figure
%     [cfolder, '\' subfolder, 'fig\', fname, '.png'])
% 
% end
% if openModel
%     open(fname); % Don't need to open to run
% end
% 
% if simulateModel
% %% Simulate System
% simConfig.SaveState = 'on';
% simOut = sim(fname, simConfig);
% 
% % Sim Data
% Xout = simOut.xout{1}.Values.Data; %Only works by grabbing states of first block (LTI_sys)
% end
% 
% if plotResults
% %% Plot Results
% fig = figure;
% plot(Xout(:,1))
% hold on
% plot(Xout(:,2))
% legend('X_1', 'X_2')
% title('Saturated Double Integration Response while stabalized')
% saveas(fig, [cfolder, '\',subfolder, 'fig\', 'StateResponse.png'])
% end

