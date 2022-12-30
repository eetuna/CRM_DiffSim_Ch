path_to_input = '3D_dynamic_response_data/input_play_files/circle_50_new.txt';
fileID = fopen(path_to_input,'r');
input_data = textscan(fileID,'%d %f %f', 'Delimiter',',');
fclose(fileID);
load('3D_dynamic_response_data/output_trajectories/circle50.mat');

circle = circle50;
Ts = 0.050; %s
Ts_camera = 0.0167; %s


index = input_data(1,1);
input_size = size(cell2mat(index), 1) / 3;

raw_currents_ = cell2mat(input_data(1,2));
raw_currents = reshape(raw_currents_, [3 input_size]);
raw_currents = 0.001*raw_currents(:,3:end); % turn it into miliamp

zero_T = floor(2/Ts);
step_T = floor(3/Ts);
current_zeros = zeros(3, zero_T);
current_x = repmat([0;0;0.2], [1 step_T]);

raw_currents = [current_zeros, current_x, raw_currents];
current_size = size(raw_currents, 2);


time_ = cell2mat(input_data(1,3));
time_1 = Ts * ones(1, zero_T + step_T);
time_input = time_(9:3:end)' * 0.001;
time_input = [time_1, time_input];
% time_input = cumsum(time_input);
total_time_input = sum(time_input);
%%outputs
cutoff_i = 1;
% cutoff_o = cutoff_i + ceil(4/ Ts_camera);
coil_position_mat = circle(cutoff_i : end,6:8) - circle(cutoff_i : end,3:5) ;
coil_position_mat = coil_position_mat' * 10^3;

normal_mat = circle(cutoff_i : end,9:11)';


time_output = circle(cutoff_i : end,1);
data_length = length(time_output);
time_output = time_output - ones(data_length, 1) *  time_output(1,1);
% 
% figure(11)
% subplot(3,2,1);
% t = Ts_camera * [1:1:length_];
% plot(t, tip_position_mat(1,1:length_), 'r-');
% subplot(3,2,3);
% t = Ts_camera * [1:1:length_];
% plot(t, tip_position_mat(2,1:length_), 'r-');
% subplot(3,2,5);
% t = Ts_camera * [1:1:length_];
% plot(t, tip_position_mat(3,1:length_), 'r-');

% subplot(3,2,2);
% plot(time_input, currents(1,1:current_size), 'k-');
% subplot(3,2,4);
% plot(time_input, currents(2,1:current_size), 'k-');
% subplot(3,2,6);
% plot(time_input, currents(3,1:current_size), 'k-');

%%output interpolation
timer = 0.0;
 
currents = zeros(3, data_length);

%interpolate inputs
for i = 1: data_length

    ind_ = floor(timer / Ts)+1; % ind_ starts from 1 not 0

    currents(:,i) = raw_currents(:, ind_);

    timer = timer + Ts_camera;
    if timer >= total_time_input
        break;
    end
end


% 
% figure(12)
% plot3(tip_position_mat(1,:), tip_position_mat(2,:), tip_position_mat(3,:), 'r-');
% hold on;
% plot3(out_put_trajectory(1,:), out_put_trajectory(2,:), out_put_trajectory(3,:), 'k*');
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');
% 
% figure(13)
% subplot(3,2,1);
% x = 1: data_length;
% plot(x, currents(1,1:data_length), 'r-');
% subplot(3,2,2);
% x = 1: data_length;
% plot(x, currents(2,1:data_length), 'r-');
% subplot(3,2,3);
% x = 1: data_length;
% plot(x, currents(3,1:data_length), 'r-');

w_currents = reshape(currents, [3*data_length 1 ]);
fileID = fopen('exp.txt','w');
formatSpec = '%f, %f, %f\n';

fprintf(fileID,formatSpec,w_currents);
