%% get estimated position and velocity
clear

Define_Constants;

%% First step
%a)
L_b = 51.509254461281560 * deg_to_rad; %latitude
lambda_b = -0.161045484887729 * deg_to_rad; %longitude
h_b = 38.825820931000635; %height
v_eb_n = [0.02;0;0]; %velocity
[r_eb_e,v_eb_e] = pv_NED_to_ECEF(L_b,lambda_b,h_b,v_eb_n);

%b)
GPSdata = load('Pseudo_ranges.csv');
[m,n] = size(GPSdata);
ranges = GPSdata(1,2:n);
t = 0;
sat_r_es_e = zeros(n-1,3);
sat_v_es_e = zeros(n-1,3);
for i = 1:n-1
[sat_r_es_e(i,:),sat_v_es_e(i,:)] = Satellite_position_and_velocity(t,ranges(:,i));
end

%c)
r_aj = zeros(n-1,1);
C_e_i = eye(3,3);
C_e_i = repmat(C_e_i,1,1,n-1);
for i = 1:n-1
    r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e)'*...
        (C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e));
    C_e_i(:,:,i) = [1,omega_ie*r_aj(i,1)/c,0;
                   -omega_ie*r_aj(i,1)/c,1,0;
                   0,0,1];
	r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e)'*...
        (C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e));
end

%Task 4 and Task 1 difference
%d)
omega_ie_e = [0,-omega_ie,0;
              omega_ie,0,0;
              0,0,0];
range_rates = load('Pseudo_range_rates.csv');

u_aj_e = zeros(n-1,3);
for i = 1:n-1
    u_aj_e(i,:) = (sat_r_es_e(i,:)' - r_eb_e)/r_aj(i,1);
end
r_aj_dot = zeros(n-1,1);
for i =1:n-1
    r_aj_dot(i,1) = u_aj_e(i,:) * (C_e_i(:,:,i) * (sat_v_es_e(i,:)' + omega_ie_e *...
        sat_r_es_e(i,:)') - (v_eb_e + omega_ie_e * r_eb_e));
end
%e)
p_a_j = GPSdata(2,2:n)'; %pseudo_range_measurement
p_a_j_dot = range_rates(2,2:n)'; %pseudo_range_rate_measurement
p_c_a = 10008.80095461663 * ones(n-1,1); %predicted receiver clock offset 
p_c_a_dot = 99.998337076036860 * ones(n-1,1); %predicted receiver clock drift

x = [r_eb_e; p_c_a(1,1)]; %pseudo_range_measurement
x_dot = [v_eb_e;p_c_a_dot(1,1)];

delta_z = p_a_j-r_aj - p_c_a; %measurement innovation vector
delta_z_dot = p_a_j_dot - r_aj_dot - p_c_a_dot;

H_g_e = [-u_aj_e ones(n-1,1)]; %measurement matrix

%f)
new_x = x + (H_g_e'*H_g_e)\H_g_e'*delta_z;
new_x_dot = x_dot + (H_g_e'*H_g_e)\H_g_e'*delta_z_dot;

%g)
[L_b_est,lambda_b_est,h_b_est,velocity] = pv_ECEF_to_NED(new_x(1:3,:),new_x_dot(1:3,:));

lattidue = L_b_est*rad_to_deg
longitude = lambda_b_est*rad_to_deg
height = h_b_est

%% Second step
L_b_est2 = zeros(11,1);
lambda_b_est2 = zeros(11,1);
h_b_est2 = zeros(11,1);
L_b_est2(1) = lattidue;
lambda_b_est2(1) = longitude;
h_b_est2(1) = height;

k = 1;
for t = GPSdata(3:852,1)'
%     %a)
%     L_b = -33.821075 * deg_to_rad; %latitude
%     lambda_b = 151.188496 * deg_to_rad; %longitude
%     h_b = 120; %height
%     v_eb_n = 0; %velocity
%     [r_eb_e,v_eb_e] = pv_NED_to_ECEF(L_b,lambda_b,h_b,v_eb_n);
    
    L_b = L_b_est2(k) * deg_to_rad; %latitude
    lambda_b = lambda_b_est2(k) * deg_to_rad; %longitude
    h_b = h_b_est2(k); %height
    v_eb_n = [0;0;0]; %velocity
    [r_eb_e,v_eb_e] = pv_NED_to_ECEF(L_b,lambda_b,h_b,v_eb_n);

    %b)
    ranges = GPSdata(1,2:n);
    sat_r_es_e = zeros(n-1,3);
    sat_v_es_e = zeros(n-1,3);
    for i = 1:n-1
    [sat_r_es_e(i,:),sat_v_es_e(i,:)] = Satellite_position_and_velocity(t,ranges(:,i));
    end

    %c)
    r_aj = zeros(n-1,1);
    C_e_i = eye(3,3);
    C_e_i = repmat(C_e_i,1,1,n-1);
    for i = 1:n-1
        r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e)'*...
            (C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e));
        C_e_i(:,:,i) = [1,omega_ie*r_aj(i,1)/c,0;
                       -omega_ie*r_aj(i,1)/c,1,0;
                       0,0,1];
        r_aj(i,1) = sqrt((C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e)'*...
            (C_e_i(:,:,i)*sat_r_es_e(i,:)' - r_eb_e));
    end

    %d)
    omega_ie_e = [0,-omega_ie,0;
              omega_ie,0,0;
              0,0,0];
    range_rates = load('Pseudo_range_rates.csv');
    u_aj_e = zeros(n-1,3);
    for i = 1:n-1
        u_aj_e(i,:) = (sat_r_es_e(i,:)' - r_eb_e)/r_aj(i,1);
    end
    for i =1:n-1
    r_aj_dot(i,1) = u_aj_e(i,:) * (C_e_i(:,:,i) * (sat_v_es_e(i,:)' + omega_ie_e *...
        sat_r_es_e(i,:)') - (v_eb_e + omega_ie_e * r_eb_e));
    end
    %e)
    p_a_j = GPSdata(k+2,2:n)'; %pseudo_range_measurement
    p_a_j_dot = range_rates(k+2,2:n)'; %pseudo_range_rate_measurement
    p_c_a = new_x(4,:) * ones(n-1,1); %predicted receiver clock offset 
    p_c_a_dot = new_x_dot(4,:) * ones(n-1,1);

    x = [r_eb_e; p_c_a(1,1)]; %pseudo_range_measurement
    x_dot = [v_eb_e;p_c_a_dot(1,1)];
    
    delta_z = p_a_j-r_aj - p_c_a; %measurement innovation vector
    delta_z_dot = p_a_j_dot - r_aj_dot - p_c_a_dot;
    
    H_g_e = [-u_aj_e ones(n-1,1)]; %measurement matrix

    %f)
    new_x = x + (H_g_e'*H_g_e)\H_g_e'*delta_z;
    new_x_dot = x_dot + (H_g_e'*H_g_e)\H_g_e'*delta_z_dot;
    
    %g)
    [L_b_est(k+1),lambda_b_est(k+1),h_b_est(k+1),velocity(:,k+1)] = pv_ECEF_to_NED(new_x(1:3,:),new_x_dot(1:3,:));

    L_b_est2(k+1) = L_b_est(k+1)*rad_to_deg;
    lambda_b_est2(k+1) = lambda_b_est(k+1)*rad_to_deg;
    h_b_est2(k+1) = h_b_est(k+1);

    k = k+1;
end    
L_b_est2(1) = lattidue;
lambda_b_est2(1) = longitude;
h_b_est2(1) = height;

plot(lambda_b_est2',L_b_est2');
csvwrite('GNSS_least_square_output.csv',[GPSdata(2:852,1),[L_b_est2,lambda_b_est2,h_b_est2,velocity']]);