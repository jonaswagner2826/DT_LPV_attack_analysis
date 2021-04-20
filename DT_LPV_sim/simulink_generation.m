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
generateModel = true;
openModel = true;
simulateModel = false;
plotResults = false;

% Name of the simulink model
[cfolder,~,~] = fileparts(mfilename('fullpath'));
subfolder = ''; %include / at end of subfolder
fname = 'DT_LPV_sim';


%% System Definitions (Simple double integrator system stabilization)
% Saturation Block
a = -1;
b = 1;

% Linearized System
A = [0, 1;
     0, 0];
B = [0;
     1];
C = eye(2);
D = [0;
     0];
lti_sys = ss(A,B,C,D);
x0 = [1;
     -1];


% Controller Gain
K = place(A,B, [-1+j,-1-j]);


%--------------------------------------------------------------------
% DT_LPV Model Definition
% (https://www.sciencedirect.com/science/article/pii/S2405896317313459)

% DT Settings
ds = 1; % Time step per k
k_max = 1000;
k = [0:ds:k_max*ds]; % Simulation Time

% Polyotopic Matrix definitions
m = 4;
A(:,:,1) = a * [-0.80, 0.25; 0.25,-0.30]; B(:,1) = [ 1.90; 0.00];
A(:,:,2) = a * [ 0.30, 0.70; 0.70, 0.00]; B(:,2) = [-1.00; 1.50];
A(:,:,3) = a * [-0.30, 0.65; 0.55, 0.10]; B(:,3) = [ 0.30;-2.00];
A(:,:,4) = a * [ 0.55,-0.20;-0.40,-0.30]; B(:,4) = [-0.60; 0.00];

A0 = [0,0;0,0]; B0 = [0;0];
A1 = [-0.80, 0.25; 0.25,-0.30]; B1 = [ 1.90; 0.00];
A2 = [ 0.30, 0.70; 0.70, 0.00]; B2 = [-1.00; 1.50];
A3 = [-0.30, 0.65; 0.55, 0.10]; B3 = [ 0.30;-2.00];
A4 = [ 0.55,-0.20;-0.40,-0.30]; B4 = [-0.60; 0.00];

% Output Eq
C = [1, 0]; D = 0;

% Sys Def
s0 = ltisys(A0,B0,C,D,0);
s1 = ltisys(A1,B1,C,D,0);
s2 = ltisys(A2,B2,C,D,0);
s3 = ltisys(A3,B3,C,D,0);
s4 = ltisys(A4,B4,C,D,0);



pv = pvec('pol',[[0;0;0;0],[0.5;0.5;0;0],[0;0;0.5;0.5],[0.5;0.5;0.5;0.5]])

poly = psys(pv,[s0,s1,s2,s3,s4])




% for i = 1:m
%     sys(i) = ss(A(:,:,i),B(:,i),C,D,ds)
% end



% sys.SamplingGrid = struct








% Input
T = 10;
u = square((2*pi*k)/10);

% Scedule Parameter
syms k_sym
g(k_sym) = piecewise(k_sym <  500, [0.50; 0.30; 0.20; 0.00],...
                     k_sym >= 500, [0.35; 0.40; 0.10; 0.15]);
alpha = g(k);

% Initial Conditions
x_0 = [0; 0];
alpha_0 = [ 0.25; 0.25; 0.25; 0.25];






if generateModel
%% Simulink Creation In Code
% Simulink Settings ----------------------
% Get the current configuration
cfg = Simulink.fileGenControl('getConfig');
% Changes Code Save Location
cfg.CacheFolder = [pwd, '\', subfolder];
cfg.CodeGenFolder = [pwd, '\', subfolder];
cfg.CodeGenFolderStructure = 'TargetEnvironmentSubfolder';
% Apply new Config
Simulink.fileGenControl('setConfig', 'config', cfg, 'keepPreviousPath',true, 'createDir',true);

% Check if the file already exists and delete it if it does
if exist(fname,'file') == 4
    % If it does then check whether it's open
    if bdIsLoaded(fname)
        % If it is then close it (without saving!)
        close_system(fname,0)
    end
    % delete the file
    delete([fname,'.slx']);
end

% Create Simulink Model
new_system;

% Create Simiple Input
add_block('simulink/Sources/In1', [gcs, '/In']);

% Create Sum block
add_block('simulink/Commonly Used Blocks/Sum', [gcs, '/Sum'],...
    'inputs', '|+-');

% Connections (Simple Input output 1 to Sum Block input 1)
add_line(gcs, 'In/1', 'Sum/1'); 

% Saturation Block
add_block('simulink/Commonly Used Blocks/Saturation', [gcs, '/Saturation'], ...
    'LowerLimit','a', ... % block parameters
    'UpperLimit','b');
add_line(gcs, 'Sum/1', 'Saturation/1');

% State-Space System
add_block('cstblocks/LTI System', [gcs, '/LTI_sys'],...
    'sys','lti_sys',...
    'IC', 'x0');
add_line(gcs, 'Saturation/1', 'LTI_sys/1');

% Controller (just a feedback gain)
add_block('simulink/Commonly Used Blocks/Gain', [gcs, '/Controller'],...
    'Gain', 'K',...
    'Multiplication', 'Matrix(K*u)',...
    'Orientation', 'left');
add_line(gcs, 'LTI_sys/1', 'Controller/1');
add_line(gcs, 'Controller/1', 'Sum/2');

% Create Simple Scope/Output
add_block('simulink/Sinks/Out1', [gcs, '/Out']);
add_line(gcs, 'LTI_sys/1', 'Out/1');



%--------------------------------------------------------------------
% LPV System

% DT_LPV_sys
add_block('cstblocks/Linear Parameter Varying/Discrete Varying State Space',...
    [gcs, '/DT_LPV_sys'])
% add_param(gcs, 'DT_LPV_sys/A', 'A(:,:,1)')









% Auto Arrange
Simulink.BlockDiagram.arrangeSystem(gcs) %Auto Arrange

%% Save and Open System
save_system(gcs,[cfolder, '\', subfolder, fname]);
print(['-s', gcs], '-dpng',... % Print model to figure
    [cfolder, '\' subfolder, 'fig\', fname, '.png'])

end
if openModel
    open(fname); % Don't need to open to run
end

if simulateModel
%% Simulate System
simConfig.SaveState = 'on';
simOut = sim(fname, simConfig);

% Sim Data
Xout = simOut.xout{1}.Values.Data; %Only works by grabbing states of first block (LTI_sys)
end

if plotResults
%% Plot Results
fig = figure;
plot(Xout(:,1))
hold on
plot(Xout(:,2))
legend('X_1', 'X_2')
title('Saturated Double Integration Response while stabalized')
saveas(fig, [cfolder, '\',subfolder, 'fig\', 'StateResponse.png'])
end

