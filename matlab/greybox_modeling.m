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


load('optimized_parameter_40ms.mat');
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

% figure(2);
% pt_r = p_mat(:,1);
% pt_d = coil_position_mat(:,ind_start);
% for i = 2: test_length
%     pt_cur_r = p_mat(:,i);
%     pt_cur_d = coil_position_mat(:,i);
% 
%     pr_ = [pt_cur_r,pt_r];
%     pd_ = [pt_cur_d,pt_d];
%     plot3(pr_(1,:), pr_(2,:), pr_(3,:), 'k-');
%     hold on;
%     plot3(pd_(1,:), pd_(2,:), pd_(3,:), 'r-');
%     hold on;
%     pt_r = pt_cur_r;
%     pt_d = pt_cur_d;
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% 
% end

% plot3(p_mat(1,:), p_mat(2,:), p_mat(3,:), 'r-');
% hold on;
% plot3(coil_position_mat(1,:), coil_position_mat(2,:), coil_position_mat(3,:), 'k*');

%% greybox
Order = [3, 3, 27]; %[Ny, Nu, Nx]

pL= [-0.531599; -1.71615; 68.9192]; %[-0.524853 ;-1.83144; 68.4151];% [-0.753845; -1.6483; 68.4184];
RL =[0.999907; -0.000299857; -0.0136261;
-0.000299857; 0.999032; -0.0439888;
0.0136261; 0.0439888 ;0.998939];
vL = [0;0;0];
wL = [0;0;0];
u0 =[0.000799849479553773; -0.0001892201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
nL = [0;0;0];
mL = [0;0;0];

% damping = [10;10;0.11;0.011]; % continuum parameter 50
% 
% radius_ = [1.5875;0.9906];
% E_ = [5.9232;1.9114]; % [5.3948;2.3881];%[5.3948;2.3881]; %[4.9105; 1.6545] 
% Coil_align =[-0.0653;-3.4378];% [-1.2; -1.5]; %[-0.0871, -0.3934]
% Coil_turnarea = [1.4;1.7005;1.2914];% [1.3851;1.44;1.60];% [1.44;1.3851;1.60];
% mass_ =[9.3516e-5];% [5.7736e-5];

% 
% damping = [10;20;0.091;0.051]; % continuum parameter 50
% radius_ = [1.5875;0.9906];
% E_ = [5.9232;1.9114]; %[4.9105; 1.6545] 
% Coil_align = [-0.0653;-3.0378]; %[-0.0871, -0.3934]
% Coil_turnarea =[1.4;1.9005;1.2914];% [1.44;1.3851;1.60];
% mass_ =[5.3516e-5];% [5.7736e-5];


load('optimized_parameter.mat');
damping = nlgr_model.Parameters(1).Value;

radius_ = [1.5875;0.9906];
E_ = nlgr_model.Parameters(4).Value;
Coil_align = nlgr_model.Parameters(5).Value; %[-0.0871, -0.3934]
Coil_turnarea = nlgr_model.Parameters(6).Value;%[1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =nlgr_model.Parameters(7).Value;% [5.7736e-5];


Parameters = {damping; Ts;radius_;E_;Coil_align;Coil_turnarea;mass_};

load('init.mat');
InitialStates = [vL,wL,u0,mL,nL,pL,RL]';

% InitialStates = [vL;wL;u0;mL;nL;pL;RL];
nlgr = idnlgrey('CRMDYN_c', Order, Parameters, InitialStates, Ts);
nlgr.Parameters(1).Fixed = false;
nlgr.Parameters(2).Fixed = true;
nlgr.Parameters(3).Fixed = true;
nlgr.Parameters(4).Fixed = false;
nlgr.Parameters(5).Fixed = false;
nlgr.Parameters(6).Fixed = false;
nlgr.Parameters(7).Fixed = false;

opt = nlgreyestOptions;
% opt.Display = 'on';
% opt.SearchOptions.MaxIterations = 2;
% opt.Advanced.ErrorThreshold=0.0001;
% opt.GradientOptions.Type = 'Basic';
opt.GradientOptions.DifferencingScheme =   'Backward approximation'; 

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

% figure(4);
% pt_r = p_mat(:,1);
% pt_d = raw_y(:,ind_y);
% for i = 2: test_length
%     pt_cur_r = p_mat(:,i);
%     pt_cur_d = raw_y(:,i+ind_y);
% 
%     pr_ = [pt_cur_r,pt_r];
%     pd_ = [pt_cur_d,pt_d];
%     plot3(pr_(1,:), pr_(2,:), pr_(3,:), 'k-');
%     hold on;
%     plot3(pd_(1,:), pd_(2,:), pd_(3,:), 'r-');
% hold on;
%     pt_r = pt_cur_r;
%     pt_d = pt_cur_d;
%     pause(0.01);
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% 
% end

ind_u = 1; % 622;
u= currents(:,ind_u:ind_u+test_length-1)';
u(:,3) = -u(:,3);
u(:,2) = -u(:,2);

z = iddata(y, u, Ts);


nlgr_model = nlgreyest(z, nlgr, opt);

% save('optimized_parameter.mat',"nlgr_model");
%% testing

pL= [-0.531599; -1.71615; 68.9192]; %[-0.524853 ;-1.83144; 68.4151];% [-0.753845; -1.6483; 68.4184];
RL =[0.999907; -0.000299857; -0.0136261;
-0.000299857; 0.999032; -0.0439888;
0.0136261; 0.0439888 ;0.998939];
vL = [0;0;0];
wL = [0;0;0];
u0 =[0.000799849479553773; -0.0002292201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
nL = [0;0;0];
mL = [0;0;0];


damping = nlgr_model.Parameters(1).Value;

radius_ = [1.5875;0.9906];
E_ = nlgr_model.Parameters(4).Value;
Coil_align = nlgr_model.Parameters(5).Value; %[-0.0871, -0.3934]
Coil_turnarea = nlgr_model.Parameters(6).Value;%[1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =nlgr_model.Parameters(7).Value;% [5.7736e-5];

load('init.mat');

ind_start = 1; %228;%622; %2
test_length = 3000; %2000; %300
p_mat = zeros(3,test_length);
for i = 1:test_length

    u = currents(:,i+ind_start);
    u(3) = -u(3);
    u(2) = -u(2);

    [vL, wL, u0, mL, nL, pL, RL] = CRMDYN_TEST_mex(vL, wL, u0, mL, nL, pL, RL, u, damping, Ts, radius_, E_,...
        Coil_align, Coil_turnarea, mass_);
    p_mat(:,i) = pL;
    pL
    RL
    i

end

ind_y = 1;

raw_y = [coil_position_mat(:,ind_y:ind_y + test_length-1)];%, normal_mat(:,ind_start:ind_start + test_length-1)'];
dy_init = coil_position_mat(:,ind_y) - p_mat(:,1);
raw_y = raw_y - repmat(dy_init, [1 test_length]);

figure(7);
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

% pause;
% 
% figure(8);
% pt_r = p_mat(:,1);
% pt_d = raw_y(:,ind_y);
% for i = 2: test_length
%     pt_cur_r = p_mat(:,i);
%     pt_cur_d = raw_y(:,i+ind_y);
% 
%     pr_ = [pt_cur_r,pt_r];
%     pd_ = [pt_cur_d,pt_d];
%     plot3(pr_(1,:), pr_(2,:), pr_(3,:), 'k-');
%     hold on;
%     plot3(pd_(1,:), pd_(2,:), pd_(3,:), 'r-');
% hold on;
%     pt_r = pt_cur_r;
%     pt_d = pt_cur_d;
%     pause(0.01);
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% 
% end
% 
% 
