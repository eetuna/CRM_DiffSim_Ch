close all;
clear;


%% input output

read_input_output;


pL= [-0.531599; -1.71615; 68.9192]; %[-0.524853 ;-1.83144; 68.4151];% [-0.753845; -1.6483; 68.4184];
RL =[0.999907; -0.000299857; -0.0136261;
-0.000299857; 0.999032; -0.0439888;
0.0136261; 0.0439888 ;0.998939];% [0.999907; -0.000324193; -0.0136272;
% -0.000324193 ;0.998869 ;-0.0475512;
% 0.0136272 ;0.0475512; 0.998776]; % [0.999808; -0.000419062; -0.019573;
%-0.000419062; 0.999084 ;-0.0427968;
%0.019573; 0.0427968 ;0.998892];
vL = [0;0;0];
wL = [0;0;0];
u0 =[0.000799849479553773; -0.0002292201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
nL = [0;0;0];
mL = [0;0;0];

% damping = [0.001;0.001;0.01;0.12;0.12;0.011];%CRM parameter
% damping = [0.001;0.001;0.01;0.1;0.1;0.01]; % continuum parameter 50
% damping = [0.0005;0.0005;0.001;0.091;0.091;0.01]; % continuum parameter 50
% damping = [0.0002;0.0002;0.001;0.08;0.08;0.01]; % continuum parameter 20, v_xy nad w_xy is different from 50

% damping = [0.0003;0.0003;0.001;0.3;0.3;0.011]; % CRM parameter 20,
% damping = [0.001;0.001;0.01;0.2;0.2;0.011]; % CRM parameter 50,

% damping = [0.001;0.01;0.2;0.011]; % CRM parameter 50,
% damping = [0.0005;0.01;0.091;0.011]; % continuum parameter 50
% damping = [0.001;0.01;0.3;0.11]; % continuum parameter 50

% damping = [10;10;0.091;0.011]; % continuum parameter 50

% damping = [0.001;0.01;0.21;0.011]; % continuum parameter 50

% damping = [10;10;0.091;0.011]; % continuum parameter 50


% radius_ = [1.5875;0.9906];
% E_ = [5.9232;1.9114]; % [5.3948;2.3881];%[5.3948;2.3881]; %[4.9105; 1.6545] 
% Coil_align =[-0.0653;-3.4378];% [-1.2; -1.5]; %[-0.0871, -0.3934]
% Coil_turnarea = [1.4;1.7005;1.2914];% [1.3851;1.44;1.60];% [1.44;1.3851;1.60];
% mass_ =[9.3516e-5];% [5.7736e-5];

load('optimized_parameter.mat');
damping = nlgr_model.Parameters(1).Value;

radius_ = [1.5875;0.9906];
E_ = nlgr_model.Parameters(4).Value;
Coil_align = nlgr_model.Parameters(5).Value; %[-0.0871, -0.3934]
Coil_turnarea = nlgr_model.Parameters(6).Value;%[1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =[9.3516e-5];% [5.7736e-5];


load('output_currents.mat');
load('output_traj.mat');
currents = output_currents;
coil_position_mat = output_traj;


load('init.mat');% at 621
ind_start = 1;%
test_length = 9500;
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
ind_start = 1;
subplot(3,1,1);
plot(x,p_mat(1, 1:test_length), 'k');
hold on;
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

damping = [10;10;0.091;0.011]; % continuum parameter 50

radius_ = [1.5875;0.9906];
E_ = [5.9232;1.9114]; % [5.3948;2.3881];%[5.3948;2.3881]; %[4.9105; 1.6545] 
Coil_align =[-0.0653;-3.4378];% [-1.2; -1.5]; %[-0.0871, -0.3934]
Coil_turnarea = [1.4;1.7005;1.2914];% [1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =[9.3516e-5];% [5.7736e-5];


Parameters = {damping; Ts;radius_;E_;Coil_align;Coil_turnarea;mass_};

load('init.mat');
InitialStates = [vL,wL,u0,mL,nL,pL,RL]';

% InitialStates = [vL;wL;u0;mL;nL;pL;RL];
nlgr = idnlgrey('CRMDYN_c', Order, Parameters, InitialStates, Ts);
nlgr.Parameters(1).Fixed = false;
nlgr.Parameters(2).Fixed = true;
nlgr.Parameters(3).Fixed = true;
nlgr.Parameters(4).Fixed = true;
nlgr.Parameters(5).Fixed = false;
nlgr.Parameters(6).Fixed = false;
nlgr.Parameters(7).Fixed = true;

opt = nlgreyestOptions;
% opt.Display = 'on';
% opt.SearchOptions.MaxIterations = 2;
% opt.Advanced.ErrorThreshold=0.0001;
% opt.GradientOptions.Type = 'Basic';
opt.GradientOptions.DifferencingScheme =   'Backward approximation'; 

ind_y = 1;

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

% % figure(4);
% % pt_r = p_mat(:,1);
% % pt_d = raw_y(:,ind_y);
% % for i = 2: test_length
% %     pt_cur_r = p_mat(:,i);
% %     pt_cur_d = raw_y(:,i+ind_y);
% % 
% %     pr_ = [pt_cur_r,pt_r];
% %     pd_ = [pt_cur_d,pt_d];
% %     plot3(pr_(1,:), pr_(2,:), pr_(3,:), 'k-');
% %     hold on;
% %     plot3(pd_(1,:), pd_(2,:), pd_(3,:), 'r-');
% % hold on;
% %     pt_r = pt_cur_r;
% %     pt_d = pt_cur_d;
% %     pause(0.01);
% % 
% %     xlabel('x');
% %     ylabel('y');
% %     zlabel('z');
% % 
% % end
% 
% ind_u = 622;
% u= currents(:,ind_u:ind_u+test_length-1)';
% u(:,3) = -u(:,3);
% u(:,2) = -u(:,2);
% 
% z = iddata(y, u, Ts);
% 
% 
% nlgr_model = nlgreyest(z, nlgr, opt);
% 
% %% testing
% 
% pL= [-0.531599; -1.71615; 68.9192]; %[-0.524853 ;-1.83144; 68.4151];% [-0.753845; -1.6483; 68.4184];
% RL =[0.999907; -0.000299857; -0.0136261;
% -0.000299857; 0.999032; -0.0439888;
% 0.0136261; 0.0439888 ;0.998939];
% vL = [0;0;0];
% wL = [0;0;0];
% u0 =[0.000799849479553773; -0.0002292201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
% nL = [0;0;0];
% mL = [0;0;0];
% 
% 
% damping = nlgr_model.Parameters(1).Value;
% 
% radius_ = [1.5875;0.9906];
% E_ = nlgr_model.Parameters(4).Value;
% Coil_align = nlgr_model.Parameters(5).Value; %[-0.0871, -0.3934]
% Coil_turnarea = nlgr_model.Parameters(6).Value;%[1.3851;1.44;1.60];% [1.44;1.3851;1.60];
% mass_ =[9.3516e-5];% [5.7736e-5];
% 
% load('init.mat');
% 
% ind_start = 622; %2
% test_length = 2000; %300
% p_mat = zeros(3,test_length);
% for i = 1:test_length
% 
%     u = currents(:,i+ind_start);
%     u(3) = -u(3);
%     u(2) = -u(2);
% 
%     [vL, wL, u0, mL, nL, pL, RL] = CRMDYN_TEST_mex(vL, wL, u0, mL, nL, pL, RL, u, damping, Ts, radius_, E_,...
%         Coil_align, Coil_turnarea, mass_);
%     p_mat(:,i) = pL;
%     pL
%     RL
%     i
% 
% end
% 
% ind_y = 500;
% 
% raw_y = [coil_position_mat(:,ind_y:ind_y + test_length-1)];%, normal_mat(:,ind_start:ind_start + test_length-1)'];
% dy_init = coil_position_mat(:,ind_y) - p_mat(:,1);
% raw_y = raw_y - repmat(dy_init, [1 test_length]);
% 
% figure(7);
% x = 1:test_length;
% 
% title('Positions');
% subplot(3,1,1);
% plot(x,p_mat(1, 1:test_length), 'k');
% hold on;
% plot(x,raw_y(1,1:1 + test_length-1), 'r');
% subplot(3,1,2);
% plot(x,p_mat(2, 1:test_length), 'k');
% hold on;
% plot(x,raw_y(2,1:1 + test_length-1), 'r');
% subplot(3,1,3);
% plot(x,p_mat(3, 1:test_length), 'k');
% hold on;
% plot(x,raw_y(3,1:1 + test_length-1), 'r');
% 
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
