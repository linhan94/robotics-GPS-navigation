function [x_est,P_matrix] = Initialise_GNSS_KF
%Initialise_GNSS_KF - Initializes the GNSS EKF state estimates and error
%covariance matrix for GNSS_KF
%
% This function created 30/11/2016 by Paul Groves
%
% Outputs:
%   x_est                 Kalman filter estimates:
%     Rows 1-3            estimated ECEF user position (m)
%     Rows 4-6            estimated ECEF user velocity (m/s)
%     Row 7               estimated receiver clock offset (m) 
%     Row 8               estimated receiver clock drift (m/s)
%   P_matrix              state estimation error covariance matrix

% Begins

% Initialise state estimates
x_est = [  3977851.131659842; -11180.86980943; 4969033.742602612;...
    0.005258490397012; 0.047195803812727;  -0.017177817280509;...
    10008.80095; 99.9983]; 

% Initialise error covariance matrix
P_matrix =  zeros(8);
P_matrix(1,1) = 100;
P_matrix(2,2) = 100;
P_matrix(3,3) = 100;
P_matrix(4,4) = 0.01;
P_matrix(5,5) = 0.01;
P_matrix(6,6) = 0.01;
P_matrix(7,7) = 100;
P_matrix(8,8) = 0.01;

% Ends