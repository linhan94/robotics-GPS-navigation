clear

Define_Constants;

[~,~,direction] = Wheel_speed;

%% Initial GNSS KF
[x_est,P_matrix] = Initialise_GNSS_KF;

%% (Kalman filter Step 1) Compute the transition matrix using
propagation_interval = 0.5;
transition_matrix = [eye(3) propagation_interval*eye(3) zeros(3,1) zeros(3,1);
                     zeros(3,3) eye(3) zeros(3,1) zeros(3,1);
                     zeros(1,3) zeros(1,3) 1 propagation_interval;
                     zeros(1,3) zeros(1,3) 0 1];

%% (Step 2) Compute the system noise covariance matrix using
acceleration_PSD = 10;
clock_phase_PSD = 1;
clock_frequency_PSD = 1;
noise_covariance_matrix = [zeros(3,3) zeros(3,3) zeros(3,1) zeros(3,1);
    zeros(3,3) acceleration_PSD*propagation_interval*eye(3) zeros(3,1) zeros(3,1);
    zeros(1,3) zeros(1,3) clock_phase_PSD*propagation_interval 0;
    zeros(1,3) zeros(1,3) 0 clock_frequency_PSD*propagation_interval]; 

%% (Step 3) Use the transition matrix to propagate the state estimates:
new_x_est = transition_matrix*x_est;

%% (Step 4) Then use this to propagate the error covariance matrix:
new_P_matrix  = transition_matrix*P_matrix*transpose(transition_matrix) + noise_covariance_matrix;

%% Compute the Cartesian ECEF positions of the satellites at time 0
pseudo_range = load('Pseudo_ranges.csv');
pseudo_range_rate = load('Pseudo_range_rates.csv');
[m,n] = size(pseudo_range);
j = pseudo_range(1,2:n);
time = 0;
sat_r_es_e = zeros(n-1,3);
sat_v_es_e = zeros(n-1,3);
for i = 1:n-1
    [sat_r_es_e(i,:),sat_v_es_e(i,:)] = Satellite_position_and_velocity(time,j(:,i));
end

r_aj = zeros(n-1,1);
C_e_i = eye(3,3);
C_e_i = repmat(C_e_i,1,1,n-1);
for i = 1:n-1
    r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:))' *...
            (C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:)));
	C_e_i(:,:,i) = [1,omega_ie*r_aj(i,1)/c,0;
                    -omega_ie*r_aj(i,1)/c,1,0;
                       0,0,1];
	r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:))' *...
            (C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:)));
end

u_aj_e = zeros(n-1,3);
for i = 1:n-1
	u_aj_e(i,:) = (sat_r_es_e(i,:)' - new_x_est(1:3,:))/r_aj(i,1);
end

%% (Step 5) Compute the measurement matrix using
H_k = [-u_aj_e,zeros(n-1,3),ones(n-1,1),zeros(n-1,1);
    zeros(n-1,3),-u_aj_e,zeros(n-1,1),ones(n-1,1)];


%% (Step 6) Compute the measurement noise covariance matrix
pseudo_range_error_std = 5;
pseudo_range_rate_error_std = 0.02;
R_k = [pseudo_range_error_std^2 * eye(n-1,n-1),zeros(n-1,n-1);
    zeros(n-1,n-1),pseudo_range_rate_error_std^2 * eye(n-1,n-1)];

%% (Step 7) Compute the Kalman gain matrix
K_k = new_P_matrix * H_k' / ((H_k * new_P_matrix * H_k') + R_k);

%% Predict the ranges and range rates from the approximate user position to each satellite
omega_ie_e = [0,-omega_ie,0;
              omega_ie,0,0;
              0,0,0];
for i =1:n-1
    r_aj_dot(i,1) = u_aj_e(i,:) * (C_e_i(:,:,i) * (sat_v_es_e(i,:)' + omega_ie_e *...
        sat_r_es_e(i,:)') - (new_x_est(4:6,:) + omega_ie_e * new_x_est(1:3,:)));
end

%% (Step 8) Formulate the measurement innovation vector
p_a_j = pseudo_range(2,2:n)'; %pseudo_range_measurement
p_a_j_dot = pseudo_range_rate(2,2:n)'; %pseudo_range_rate_measurement
p_c_a = new_x_est(7,:); %propagated receiver clock offset estimate 
p_c_a_dot = new_x_est(8,:); %propagated receiver clock offset estimate

delta_z = [p_a_j - r_aj - p_c_a;
            p_a_j_dot - r_aj_dot - p_c_a_dot];%measurement innovation vector

%% (Step 9) Update the state estimates using
Updated_x_est = new_x_est + K_k * delta_z;

%% (Step 10) Update the error covariance matrix using
Updated_P_matrix = (eye(8) - K_k * H_k) * new_P_matrix;

%% Convert this Cartesian ECEF position solution to latitude, longitude and height
[L_b_est,lambda_b_est,h_b_est,velocity] = pv_ECEF_to_NED(Updated_x_est(1:3,:),Updated_x_est(4:6,:));
lattidue = L_b_est*rad_to_deg
longitude = lambda_b_est*rad_to_deg
height = h_b_est

k = 1;
%% Get result for the following time
for t = pseudo_range(3:852,1)'
    %% Initial GNSS KF
    x_est = Updated_x_est;
    P_matrix = Updated_P_matrix;

    %% (Kalman filter Step 1) Compute the transition matrix using
    if t == 34 || t == 35
        propagation_interval = 1;
    else
        propagation_interval = 0.5;
    end
    transition_matrix = [eye(3) propagation_interval*eye(3) zeros(3,1) zeros(3,1);
                         zeros(3,3) eye(3) zeros(3,1) zeros(3,1);
                         zeros(1,3) zeros(1,3) 1 propagation_interval;
                         zeros(1,3) zeros(1,3) 0 1];

    %% (Step 2) Compute the system noise covariance matrix using
    acceleration_PSD = 10;
    clock_phase_PSD = 1;
    clock_frequency_PSD = 1;
    noise_covariance_matrix = [zeros(3,3) zeros(3,3) zeros(3,1) zeros(3,1);
        zeros(3,3) acceleration_PSD*propagation_interval*eye(3) zeros(3,1) zeros(3,1);
        zeros(1,3) zeros(1,3) clock_phase_PSD*propagation_interval 0;
        zeros(1,3) zeros(1,3) 0 clock_frequency_PSD*propagation_interval]; 

    %% (Step 3) Use the transition matrix to propagate the state estimates:
    new_x_est = transition_matrix*x_est;

    %% (Step 4) Then use this to propagate the error covariance matrix:
    new_P_matrix  = transition_matrix*P_matrix*transpose(transition_matrix) + noise_covariance_matrix;

    %% Compute the Cartesian ECEF positions of the satellites at time 0
    j = pseudo_range(1,2:n);
    sat_r_es_e = zeros(n-1,3);
    sat_v_es_e = zeros(n-1,3);
    for i = 1:n-1
        [sat_r_es_e(i,:),sat_v_es_e(i,:)] = Satellite_position_and_velocity(t,j(:,i));
    end

    r_aj = zeros(n-1,1);
    C_e_i = eye(3,3);
    C_e_i = repmat(C_e_i,1,1,n-1);
    for i = 1:n-1
        r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:))' *...
                (C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:)));
        C_e_i(:,:,i) = [1,omega_ie*r_aj(i,1)/c,0;
                        -omega_ie*r_aj(i,1)/c,1,0;
                           0,0,1];
        r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:))' *...
                (C_e_i(:,:,i)*sat_r_es_e(i,:)' - new_x_est(1:3,:)));
    end

    u_aj_e = zeros(n-1,3);
    for i = 1:n-1
        u_aj_e(i,:) = (sat_r_es_e(i,:)' - new_x_est(1:3,:))/r_aj(i,1);
    end

    %% (Step 5) Compute the measurement matrix using
    H_k = [-u_aj_e,zeros(n-1,3),ones(n-1,1),zeros(n-1,1);
        zeros(n-1,3),-u_aj_e,zeros(n-1,1),ones(n-1,1)];


    %% (Step 6) Compute the measurement noise covariance matrix
    pseudo_range_error_std = 205.2;
    pseudo_range_rate_error_std = 0.02;
    R_k = [pseudo_range_error_std^2 * eye(n-1,n-1),zeros(n-1,n-1);
        zeros(n-1,n-1),pseudo_range_rate_error_std^2 * eye(n-1,n-1)];

    %% (Step 7) Compute the Kalman gain matrix
    K_k = new_P_matrix * H_k' / ((H_k * new_P_matrix * H_k') + R_k);

    %% Predict the ranges and range rates from the approximate user position to each satellite
    omega_ie_e = [0,-omega_ie,0;
                  omega_ie,0,0;
                  0,0,0];
    for i =1:n-1
        r_aj_dot(i,1) = u_aj_e(i,:) * (C_e_i(:,:,i) * (sat_v_es_e(i,:)' + omega_ie_e *...
            sat_r_es_e(i,:)') - (new_x_est(4:6,:) + omega_ie_e * new_x_est(1:3,:)));
    end

    %% (Step 8) Formulate the measurement innovation vector
    p_a_j = pseudo_range(k+2,2:n)'; %pseudo_range_measurement
    p_a_j_dot = pseudo_range_rate(k+2,2:n)'; %pseudo_range_rate_measurement
    p_c_a = new_x_est(7,:); %propagated receiver clock offset estimate 
    p_c_a_dot = new_x_est(8,:); %propagated receiver clock offset estimate

    delta_z = [p_a_j - r_aj - p_c_a;
                p_a_j_dot - r_aj_dot - p_c_a_dot];%measurement innovation vector

    %% (Step 9) Update the state estimates using
    Updated_x_est = new_x_est + K_k * delta_z;

    %% (Step 10) Update the error covariance matrix using
    Updated_P_matrix = (eye(8) - K_k * H_k) * new_P_matrix;

    %% Convert this Cartesian ECEF position solution to latitude, longitude and height
    [L_b_est(k+1),lambda_b_est(k+1),h_b_est(k+1),velocity(:,k+1)] = pv_ECEF_to_NED(Updated_x_est(1:3,:),Updated_x_est(4:6,:));
    L_b_est2(k+1) = L_b_est(k+1)*rad_to_deg;
    lambda_b_est2(k+1) = lambda_b_est(k+1)*rad_to_deg;
    h_b_est2(k+1) = h_b_est(k+1);
    
    k = k+1;
end

L_b_est2(1) = lattidue;
lambda_b_est2(1) = longitude;
h_b_est2(1) = height;

figure
plot(lambda_b_est2,L_b_est2);
title('Pure GNSS Navigation');
xlabel('longitude');
ylabel('lattidue');

dlmwrite('GNSS_KF_output.csv',[pseudo_range(2:852,1),L_b_est2',lambda_b_est2',velocity(1:2,:)',direction],'precision',15);