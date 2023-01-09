path_to_input = '3D_dynamic_response_data/input_play_files/circle_100_new.txt';
fileID = fopen(path_to_input,'r');
input_data = textscan(fileID,'%d %f %f', 'Delimiter',',');
fclose(fileID);
load('3D_dynamic_response_data/output_trajectories/circle100.mat');

circle = circle100;
Ts = 0.10; %s
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
coil_position_mat_raw = circle(cutoff_i : end,6:8) - circle(cutoff_i : end,3:5) ;
coil_position_mat_raw = coil_position_mat_raw' * 10^3;

normal_mat = circle(cutoff_i : end,9:11)';


time_output = circle(cutoff_i : end,1);
data_length = length(time_output);
time_output = [0:1:data_length-1] * Ts_camera; %time_output - ones(data_length, 1) *  time_output(1,1);

% figure(11)
% length_ = size(coil_position_mat_raw,2);
% subplot(3,2,1);
% t = [1:1:length_];
% plot(t, coil_position_mat_raw(1,1:length_), 'r-');
% subplot(3,2,3);
% t = [1:1:length_];
% plot(t, coil_position_mat_raw(2,1:length_), 'r-');
% subplot(3,2,5);
% t = [1:1:length_];
% plot(t, coil_position_mat_raw(3,1:length_), 'r-');
% 
% subplot(3,2,2);
% t = [1:1:current_size];
% plot(t, raw_currents(1,1:current_size), 'k-');
% subplot(3,2,4);
% plot(t, raw_currents(2,1:current_size), 'k-');
% subplot(3,2,6);
% plot(t, raw_currents(3,1:current_size), 'k-');

%%output interpolation

if Ts >= Ts_camera
%     timer = 0.0; 
%     currents = zeros(3, data_length);
%     %interpolate inputs
%     for i = 1: data_length
%     
%         ind_ = floor(timer / Ts)+1; % ind_ starts from 1 not 0
%     
%         currents(:,i) = raw_currents(:, ind_);
%     
%         timer = timer + Ts_camera;
%         if timer >= total_time_input
%             break;
%         end
%     end
%     coil_position_mat = coil_position_mat_raw;

    coil_position_mat = [];
    timer = 0.0;
    
    for i = 1: current_size
        ind_up = ceil(timer / Ts_camera)+1;
        ind_down = floor(timer / Ts_camera)+1;
        interpolated_pts = (timer - time_output(ind_down)) * (coil_position_mat_raw(:,ind_up)...
            - coil_position_mat_raw(:,ind_down)) / Ts_camera + coil_position_mat_raw(:,ind_down);
    
        coil_position_mat = [coil_position_mat, interpolated_pts];
        timer = timer + Ts;
        if timer >= time_output(data_length)
            break;
        end
    end
    data_length = current_size;
    currents = raw_currents(:, 1:current_size);

    figure(12)
    subplot(3,2,1);
    length_ = size(coil_position_mat,2);
    t = [1:1:length_] ;
    plot(t, coil_position_mat(1,1:length_), 'r-');
    hold on;
    length_ = size(coil_position_mat_raw,2);
    t = [1:1:length_] * Ts_camera;
    plot(t, coil_position_mat_raw(1,1:length_), 'b-');
    subplot(3,2,3);
    length_ = size(coil_position_mat,2);
    t = [1:1:length_]* Ts;
    plot(t, coil_position_mat(2,1:length_), 'r-');
    hold on;
    length_ = size(coil_position_mat_raw,2);
    t = [1:1:length_]* Ts_camera;
    plot(t, coil_position_mat_raw(2,1:length_), 'b-');
    subplot(3,2,5);
     length_ = size(coil_position_mat,2);
    t = [1:1:length_]* Ts;
    plot(t, coil_position_mat(3,1:length_), 'r-');
    hold on;
    length_ = size(coil_position_mat_raw,2);
    t = [1:1:length_]* Ts_camera;
    plot(t, coil_position_mat_raw(3,1:length_), 'b-');
    
    subplot(3,2,2);
    t = [1:1:data_length];
    plot(t, currents(1,1:data_length), 'k-');
    subplot(3,2,4);
    plot(t, currents(2,1:data_length), 'k-');
    subplot(3,2,6);
    plot(t, currents(3,1:data_length), 'k-');



%     figure(12)
%     subplot(3,1,1)
%     t = [1:1:data_length]*Ts_camera;
%     plot(t, currents(1,1:data_length), 'k-');
%     hold on;
%     t = [1:1:current_size];
%     plot(t, raw_currents(1,1:current_size), 'r-'); 
%     subplot(3,1,2);
%     t = [1:1:data_length];
% 
%     plot(t, currents(2,1:data_length), 'k-');
%     hold on;
%     t = [1:1:current_size];
%     plot(t, raw_currents(2,1:current_size), 'r-'); 
%     subplot(3,1,3);
%     t = [1:1:data_length]* Ts;
%     plot(t, currents(3,1:data_length), 'k-');
%     hold on;
%     t = [1:1:current_size];
%     plot(t, raw_currents(3,1:current_size), 'r-'); 
end

if Ts < Ts_camera

    coil_position_mat = [];
    timer = 0.0;
    
    for i = 1: data_length
        ind_up = ceil(timer / Ts_camera)+1;
        ind_down = floor(timer / Ts_camera)+1;
        interpolated_pts = (timer - time_output(ind_down)) * (coil_position_mat_raw(:,ind_up)...
            - coil_position_mat_raw(:,ind_down)) / Ts_camera + coil_position_mat_raw(:,ind_down);
    
        coil_position_mat = [coil_position_mat, interpolated_pts];
        timer = timer + Ts;
        if timer >= time_output(data_length)
            break;
        end
    end
    
    data_length = size(coil_position_mat, 2);
    currents = raw_currents(:, 1:data_length);

    figure(12)
    length_ = size(coil_position_mat,2);
    subplot(3,2,1);
    t = [1:1:length_];
    plot(t, coil_position_mat(1,1:length_), 'r-');
    hold on;
    length_ = size(coil_position_mat_raw,2);
    t = [1:1:length_] * Ts_camera;
    plot(t, coil_position_mat_raw(1,1:length_), 'b-');
    subplot(3,2,3);
    t = [1:1:length_]* Ts;
    plot(t, coil_position_mat(2,1:length_), 'r-');
    hold on;
    length_ = size(coil_position_mat_raw,2);
    t = [1:1:length_]* Ts_camera;
    plot(t, coil_position_mat_raw(2,1:length_), 'b-');
    subplot(3,2,5);
    t = [1:1:length_]* Ts;
    plot(t, coil_position_mat(3,1:length_), 'r-');
    hold on;
    length_ = size(coil_position_mat_raw,2);
    t = [1:1:length_]* Ts_camera;
    plot(t, coil_position_mat_raw(3,1:length_), 'b-');
    
    subplot(3,2,2);
    t = [1:1:data_length];
    plot(t, currents(1,1:data_length), 'k-');
    subplot(3,2,4);
    plot(t, currents(2,1:data_length), 'k-');
    subplot(3,2,6);
    plot(t, currents(3,1:data_length), 'k-');


end