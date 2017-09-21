 function final_direction = direction
% %DIRECTION Summary of this function goes here
% %   Detailed explanation goes here
    Define_Constants;
    DRdata = load('Dead_reckoning.csv');

    time = DRdata(:,1);
    gyroData = DRdata(:,6);
    magneticheadingData = DRdata(:,7);
    trueheadingData = DRdata(:,7) - ones(size(magneticheadingData));
    
    final_direction = zeros(851,1);
    final_direction(1,:) = deg_to_rad*trueheadingData(1);
    %low pass filter
%     filtered_gyroData = zeros(851,1);
%     for i = 1:851
%         if abs(gyroData(i,:)) > 0.04
%             filtered_gyroData(i) = gyroData(i,:);
%         end
%     end
    
    k = 0;
    for i = 2:851
        final_direction(i,:) = final_direction(i-1,:) + k * gyroData(i-1)...
        * (time(i)-time(i-1)) + deg_to_rad*(1-k) * (trueheadingData(i)-trueheadingData(i-1));
    end
 end
% 
