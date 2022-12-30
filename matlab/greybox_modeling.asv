close all;
clear;


%% input output

read_input_output;

%% preparation for inputs
vL = [0;0;0];
wL = [0;0;0];


u0 =[0.0007399849479553773; -0.0001892201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
n0 = [0;0;0];
uL = [0.000148626102272827; 0.00094448815853795; 0 ];
nL = [0;0;0];

damping_tubing = [0.5;0.1;0.1];%CRM parameter
% damping_coil = [0.0005;0.01;0.01;0.0001]; % continuum parameter 50
damping_coil = [0.0005;0.001;0.01;0.0001]; % continuum parameter 50



radius_ = [1.5875;0.9906];

%parameters that worked...

% E_ = [5.2232;2.1114]; 
% Coil_align =[-0.0653;-3.0378];
% Coil_turnarea = [1.4;1.6005;1.1914];%
% mass_ =[9.3516e-5];% [5.7736e-5];

E_ = [5.10232;2.1114]; % [5.3948;2.3881];%[5.3948;2.3881]; %[4.9105; 1.6545] 
Coil_align =[-0.0753;-2.90378];% [-1.2; -1.5]; %[-0.0871, -0.3934]
Coil_turnarea = [1.4;1.5005;1.1914];% [1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =[9.3516e-5];% [5.7736e-5];

IntegrationStepSize = 0.2;
SegmentLengths = [19.85;  18.3; 59.40];
ustarlist  = [0.000148626102272827; 0.00094448815853795; 0; 0.000799849479553773; -0.0001892201697658481; 0];

NUM_SEGMENTS  = 3;
NUM_FLEX_SEG = 2;
SegEnds = zeros(NUM_SEGMENTS+1, 1); %[NUM_SEGMENTS+1];
SegEnds(1) = 0.0;
for i = 2:NUM_SEGMENTS+1
    SegEnds(i) = SegEnds(i-1) + SegmentLengths(NUM_SEGMENTS - i+2);
end

SegSteps = zeros(NUM_FLEX_SEG, 1);
h0 = zeros(NUM_FLEX_SEG, 1);
for i =1: NUM_FLEX_SEG
    temp = ceil( (SegEnds(2*i)-SegEnds(2*i-1) ) / IntegrationStepSize );
    SegSteps(i) = int64(temp);
    h0(i) =( SegEnds(2*i)-SegEnds(2*i-1) )/(SegSteps(i)*1.0);

end

length_1 = ( SegSteps(1) + 1 ) * 3 ;
length_2 = ( SegSteps(2) + 1 ) * 3 ;
fseg_lengths = [length_1; length_2];

u_history = zeros(length_1 + length_2, 1);
for i = 1: SegSteps(1) +1
    u_history(3 *(i-1)+1) = ustarlist(4);
    u_history(3 *(i-1)+2) = ustarlist(5);
    u_history(3 *(i-1)+3) = ustarlist(6);
end
for i = SegSteps(1) + 1 + 1 : SegSteps(1) + 1 + SegSteps(2) + 1
    u_history(3 *(i-1)+1) = ustarlist(1);
    u_history(3 *(i-1)+2) = ustarlist(2);
    u_history(3 *(i-1)+3) = ustarlist(3);
end

v_history = zeros(length_1 + length_2, 1);
for i = 1: SegSteps(1) +1
    v_history(3 *(i-1)+1) = 0.0;
    v_history(3 *(i-1)+2) = 0.0;
    v_history(3 *(i-1)+3) = 0.0;
end
for i = SegSteps(1) + 1 + 1 : SegSteps(1) + 1 + SegSteps(2) + 1
    v_history(3 *(i-1)+1) = 0.0;
    v_history(3 *(i-1)+2) = 0.0;
    v_history(3 *(i-1)+3) = 0.0;
end

w_history = zeros(length_1 + length_2, 1);
for i = 1: SegSteps(1) +1
    w_history(3 *(i-1)+1) = 0.0;
    w_history(3 *(i-1)+2) = 0.0;
    w_history(3 *(i-1)+3) = 0.0;
end
for i = SegSteps(1) + 1 + 1 : SegSteps(1) + 1 + SegSteps(2) + 1
    w_history(3 *(i-1)+1) = 0.0;
    w_history(3 *(i-1)+2) = 0.0;
    w_history(3 *(i-1)+3) = 0.0;
end

load('init.mat');

ind_start = 621;%2;
test_length = 1200;
p_mat = [];

for i = ind_start:ind_start + test_length-1

    u = currents(:,i);
    u(3) = -u(3);  
    u(2) = -u(2);
%     temp_ = u(1);
%     u(1) = u(2);
%     u(2) = temp_;

[vL, wL, u0, n0, uL, nL, pL, RL, u_history, v_history, w_history] = CRMDYN_TEST_mex(vL, wL, u0, n0, uL, nL, u, ...
    damping_tubing, damping_coil, Ts, radius_, E_, Coil_align, Coil_turnarea, SegmentLengths,...
    mass_, IntegrationStepSize, SegSteps, u_history, v_history, w_history);


    
    p_mat = [p_mat,pL'];

    vL
    wL
    pL
    RL


    i

end
figure(1);
ind_start=500;
x = 1:test_length;

title('Positions');

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

% save('init.mat','vL', 'wL', 'u0', 'n0', 'uL', 'nL', 'pL', 'RL', 'h0', 'u_history', 'v_history', 'w_history');

% 
% 
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
%     pause(0.01);
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% 
% end
% plot3(p_mat(1,:), p_mat(2,:), p_mat(3,:), 'r-');
% hold on;
% plot3(coil_position_mat(1,:), coil_position_mat(2,:), coil_position_mat(3,:), 'k*');
% 
% 

%% greybox

vL = [0;0;0];
wL = [0;0;0];

SegEnds = zeros(NUM_SEGMENTS+1, 1); %[NUM_SEGMENTS+1];
SegEnds(1) = 0.0;
for i = 2:NUM_SEGMENTS+1
    SegEnds(i) = SegEnds(i-1) + SegmentLengths(NUM_SEGMENTS - i+2);
end

SegSteps = zeros(NUM_FLEX_SEG, 1);
h0 = zeros(NUM_FLEX_SEG, 1);
for i =1: NUM_FLEX_SEG
    temp = ceil( (SegEnds(2*i)-SegEnds(2*i-1) ) / IntegrationStepSize );
    SegSteps(i) = int64(temp);
end


u0 =[0.0007399849479553773; -0.0001892201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
n0 = [0;0;0];
uL = [0.000148626102272827; 0.00094448815853795; 0 ];
nL = [0;0;0];

damping_tubing = [0.1;0.1;0.1];%CRM parameter
% damping_coil = [0.0005;0.01;0.01;0.0001]; % continuum parameter 50
damping_coil = [0.0005;0.0001;0.1;0.0001]; % continuum parameter 50


radius_ = [1.5875;0.9906];
E_ = [5.10232;2.1114]; % [5.3948;2.3881];%[5.3948;2.3881]; %[4.9105; 1.6545] 
Coil_align =[-0.0753;-2.90378];% [-1.2; -1.5]; %[-0.0871, -0.3934]
Coil_turnarea = [1.4;1.5005;1.1914];% [1.3851;1.44;1.60];% [1.44;1.3851;1.60];
mass_ =[9.3516e-5];% [5.7736e-5];
IntegrationStepSize = 0.2;
SegmentLengths = [19.85;  18.3; 59.40];

Nx = 3*7 +9 + 3*(length_1 + length_2);
Order = [3, 3, Nx]; %[Ny, Nu, Nx]

Parameters = {SegSteps; damping_tubing; damping_coil;Ts;radius_;E_;Coil_align;Coil_turnarea;SegmentLengths;mass_;IntegrationStepSize};


pL = [0;0;0]; % not used in the calculation
RL=[1;0;0;0;1;0;0;0;1];

load('init.mat');

InitialStates = [vL,wL,u0,n0,uL,nL,pL,RL,u_history,v_history,w_history]';


nlgr = idnlgrey('CRMDYN_c', Order, Parameters, InitialStates, Ts);
nlgr.Parameters(1).Fixed = true;
nlgr.Parameters(2).Fixed = false;
nlgr.Parameters(3).Fixed = false;
nlgr.Parameters(4).Fixed = true;
nlgr.Parameters(5).Fixed = true;
nlgr.Parameters(6).Fixed = false;
nlgr.Parameters(7).Fixed = false;
nlgr.Parameters(8).Fixed = false;
nlgr.Parameters(9).Fixed = true;
nlgr.Parameters(10).Fixed = true;
nlgr.Parameters(11).Fixed = true;

opt = nlgreyestOptions;
% opt.Display = 'on';
% opt.SearchOptions.MaxIterations = 2;
% opt.Advanced.ErrorThreshold= 10;
% opt.GradientOptions.Type = 'Basic';
opt.GradientOptions.DifferencingScheme =   'Backward approximation'; 


ind_y = 500;

raw_y = [coil_position_mat(:,ind_y:ind_y + test_length-1)];%, normal_mat(:,ind_start:ind_start + test_length-1)'];
dy_init = coil_position_mat(:,ind_y) - p_mat(:,1);
raw_y = raw_y - repmat(dy_init, [1 test_length]);
y = raw_y';
% y = p_mat';
% 
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


ind_u = 621;
u= currents(:,ind_u:ind_u + test_length-1)';
u(:,3) = -u(:,3);
u(:,2) = -u(:,2);
% temp_ = u(:,1);
% u(:,1) = u(:,2);
% u(:,2) = temp_;

z = iddata(y, u, Ts);

nlgr_model = nlgreyest(z, nlgr, opt);

% %% testing
% vL = [0;0;0];
% wL = [0;0;0];
% 
% 
% u0 =[0.0007399849479553773; -0.0001892201697658481; 0];% [0.000719849479553773; -0.0003292201697658481; 0];
% n0 = [0;0;0];
% uL = [0.000148626102272827; 0.00094448815853795; 0 ];
% nL = [0;0;0];
% 
% damping_tubing = nlgr_model.Parameters(2).Value
% damping_coil =  nlgr_model.Parameters(3).Value
% 
% 
% radius_ = [1.5875;0.9906];
% E_ = nlgr_model.Parameters(6).Value
% Coil_align = nlgr_model.Parameters(7).Value
% Coil_turnarea = nlgr_model.Parameters(8).Value
% mass_ =[9.3516e-5];% [5.7736e-5];
% 
% IntegrationStepSize = 0.2;
% SegmentLengths = [19.85;  18.3; 59.40];
% ustarlist  = [0.000148626102272827; 0.00094448815853795; 0; 0.000799849479553773; -0.0001892201697658481; 0];
% 
% NUM_SEGMENTS  = 3;
% NUM_FLEX_SEG = 2;
% SegEnds = zeros(NUM_SEGMENTS+1, 1); %[NUM_SEGMENTS+1];
% SegEnds(1) = 0.0;
% for i = 2:NUM_SEGMENTS+1
%     SegEnds(i) = SegEnds(i-1) + SegmentLengths(NUM_SEGMENTS - i+2);
% end
% 
% SegSteps = zeros(NUM_FLEX_SEG, 1);
% h0 = zeros(NUM_FLEX_SEG, 1);
% for i =1: NUM_FLEX_SEG
%     temp = ceil( (SegEnds(2*i)-SegEnds(2*i-1) ) / IntegrationStepSize );
%     SegSteps(i) = int64(temp);
%     h0(i) =( SegEnds(2*i)-SegEnds(2*i-1) )/(SegSteps(i)*1.0);
% 
% end
% 
% length_1 = ( SegSteps(1) + 1 ) * 3 ;
% length_2 = ( SegSteps(2) + 1 ) * 3 ;
% fseg_lengths = [length_1; length_2];
% 
% u_history = zeros(length_1 + length_2, 1);
% for i = 1: SegSteps(1) +1
%     u_history(3 *(i-1)+1) = ustarlist(4);
%     u_history(3 *(i-1)+2) = ustarlist(5);
%     u_history(3 *(i-1)+3) = ustarlist(6);
% end
% for i = SegSteps(1) + 1 + 1 : SegSteps(1) + 1 + SegSteps(2) + 1
%     u_history(3 *(i-1)+1) = ustarlist(1);
%     u_history(3 *(i-1)+2) = ustarlist(2);
%     u_history(3 *(i-1)+3) = ustarlist(3);
% end
% 
% v_history = zeros(length_1 + length_2, 1);
% for i = 1: SegSteps(1) +1
%     v_history(3 *(i-1)+1) = 0.0;
%     v_history(3 *(i-1)+2) = 0.0;
%     v_history(3 *(i-1)+3) = 0.0;
% end
% for i = SegSteps(1) + 1 + 1 : SegSteps(1) + 1 + SegSteps(2) + 1
%     v_history(3 *(i-1)+1) = 0.0;
%     v_history(3 *(i-1)+2) = 0.0;
%     v_history(3 *(i-1)+3) = 0.0;
% end
% 
% w_history = zeros(length_1 + length_2, 1);
% for i = 1: SegSteps(1) +1
%     w_history(3 *(i-1)+1) = 0.0;
%     w_history(3 *(i-1)+2) = 0.0;
%     w_history(3 *(i-1)+3) = 0.0;
% end
% for i = SegSteps(1) + 1 + 1 : SegSteps(1) + 1 + SegSteps(2) + 1
%     w_history(3 *(i-1)+1) = 0.0;
%     w_history(3 *(i-1)+2) = 0.0;
%     w_history(3 *(i-1)+3) = 0.0;
% end
% 
% load('init.mat');
% 
% ind_start = 621;%2;
% test_length = 2000;
% p_mat = [];
% 
% for i = ind_start:ind_start + test_length-1
% 
%     u = currents(:,i);
%     u(3) = -u(3);  
%     u(2) = -u(2);
% 
% [vL, wL, u0, n0, uL, nL, pL, RL, h0, u_history, v_history, w_history] = CRMDYN_TEST_mex(vL, wL, u0, n0, uL, nL, u, ...
%     damping_tubing, damping_coil, Ts, radius_, E_, Coil_align, Coil_turnarea, SegmentLengths,...
%     mass_, IntegrationStepSize, SegSteps, h0, u_history, v_history, w_history);
% 
% 
%     
%     p_mat = [p_mat,pL'];
% 
%     vL
%     wL
%     pL
%     RL
% 
% 
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
% % pause;
% % 
% % figure(8);
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