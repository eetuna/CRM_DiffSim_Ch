close all;
clear;


%% input output

read_input_output;


pL= [-0.531599; -1.71615; 68.9192]; %[-0.524853 ;-1.83144; 68.4151];% [-0.753845; -1.6483; 68.4184];
RL =[0.999907; -0.000299857; -0.0136261;
-0.000299857; 0.999032; -0.0439888;
0.0136261; 0.0439888 ;0.998939];

vL = [0;0;0];
wL = [0;0;0];
u0 =[0.000799849479553773; -0.0002292201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
nL = [0;0;0];
mL = [0;0;0];

% damping = [10;10;0.091;0.011]; % continuum parameter 50
% radius_ = [1.5875;0.9906];
% E_ = [5.9232;1.9114]; % [5.3948;2.3881];%[5.3948;2.3881]; %[4.9105; 1.6545] 
% Coil_align =[-0.0653;-3.4378];% [-1.2; -1.5]; %[-0.0871, -0.3934]
% Coil_turnarea = [1.4;1.7005;1.2914];% [1.3851;1.44;1.60];% [1.44;1.3851;1.60];
% mass_ =[9.3516e-5];% [5.7736e-5];
% 
% continuum parameter 50
% damping = [10;20;0.091;0.051]; 
% radius_ = [1.5875;0.9906];
% E_ = [5.9232;1.9114]; %[4.9105; 1.6545] 
% Coil_align = [-0.0653;-3.0378]; %[-0.0871, -0.3934]
% Coil_turnarea =[1.4;1.9005;1.2914];% [1.44;1.3851;1.60];
% mass_ =[5.3516e-5];% [5.7736e-5];


% % continuum parameter 100
% damping = [10;20;0.1;0.051]; 
% radius_ = [1.5875;0.9906];
% E_ = [5.9232;1.9114]; %[4.9105; 1.6545] 
% Coil_align = [-0.0653;-3.0378]; %[-0.0871, -0.3934]
% Coil_turnarea =[1.4;1.9005;1.2914];% [1.44;1.3851;1.60];
% mass_ =[5.3516e-5];% [5.7736e-5];


load(['optimized_parameter_40ms.mat']);
damping = nlgr_model.Parameters(1).Value;

radius_ = [1.5875;0.9906];
E_ = nlgr_model.Parameters(4).Value;
Coil_align = nlgr_model.Parameters(5).Value; %[-0.0871, -0.3934]
Coil_turnarea = nlgr_model.Parameters(6).Value;%[1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =nlgr_model.Parameters(7).Value;% [5.7736e-5];

%test with varying freq data
load('output_currents.mat');
load('output_traj.mat');
currents = output_currents;
coil_position_mat = output_traj;
Ts = 0.04;

load('init.mat');% at 621
ind_start = 1;% 228 622;
test_length = 3000; %621
p_mat = [];
for i = ind_start:ind_start+test_length-1

    u = currents(:,i);
    u(3) = -u(3);
    u(2) = -u(2);
%     temp_ = u(1);
%     u(1) = u(2);
%     u(2) = temp_;

    [vL, wL, u0, mL, nL, pL, RL] = CRMDYN_TEST_mex(vL, wL, u0, mL, nL, pL, RL, u, damping, Ts, radius_, E_,...
        Coil_align, Coil_turnarea, mass_);
       
    p_mat = [p_mat,pL'];
    pL
    RL
    i

end

% save('init.mat','vL','wL','u0','mL','nL','pL','RL');

figure(1);
x = 1:test_length;

title('Positions');
ind_start = 1; %142; % 171;500;
subplot(3,1,1);
x = [1:test_length];
plot(x,p_mat(1, 1:test_length), 'k');
hold on;
x = [1:test_length] ;
plot(x,coil_position_mat(1,ind_start:ind_start + test_length-1), 'r');
subplot(3,1,2);
plot(x,p_mat(2, 1:test_length), 'k');
hold on;
plot(x,coil_position_mat(2,ind_start:ind_start + test_length-1), 'r');
subplot(3,1,3);
plot(x,p_mat(3, 1:test_length), 'k');
hold on;
plot(x,coil_position_mat(3,ind_start:ind_start + test_length-1), 'r');



ind_y = 1; %142; %171;

raw_y = [coil_position_mat(:,ind_y:ind_y + test_length-1)];%, normal_mat(:,ind_start:ind_start + test_length-1)'];
dy_init = coil_position_mat(:,ind_y) - p_mat(:,1);
raw_y = raw_y - repmat(dy_init, [1 test_length]);
y = raw_y';
% y = p_mat';

figure(3);
x = 1:test_length;

title('Positions');
subplot(3,1,1);
plot(x,p_mat(1, 1:test_length), 'k');
hold on;
plot(x,raw_y(1,1:1 + test_length-1), 'r');
subplot(3,1,2);
plot(x,p_mat(2, 1:test_length), 'k');
hold on;
plot(x,raw_y(2,1:1 + test_length-1), 'r');
subplot(3,1,3);
plot(x,p_mat(3, 1:test_length), 'k');
hold on;
plot(x,raw_y(3,1:1 + test_length-1), 'r');
