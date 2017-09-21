function [position,v_eb,vehicleDirection] = Wheel_speed
    vehicleDirection = direction;

    DRdata = load('Dead_reckoning.csv');
    Define_Constants
    
    l_bf_b = 0.3;
    l_br_b = -0.2;
    time = DRdata(:,1);
    wheel_speed = DRdata(:,2:5);
    [m,n] = size(wheel_speed);
    %The locations of the wheel speed sensors with respect to the GNSS antenna are as follows:
    %1) 0.3m forwards, -0.2m right;
    %2) 0.3m forwards, 0.2m right;
    %3) -0.2m forwards, -0.2m right;
    %4) -0.2m forwards, 0.2m right.

    v_efL_f = [wheel_speed(:,1)';
               zeros(1,m)];
    v_efR_f = [wheel_speed(:,2)';
               zeros(1,m)];
    v_erL_r = [wheel_speed(:,3)';
               zeros(1,m)];
    v_erR_r = [wheel_speed(:,4)';
               zeros(1,m)];

    v_er = (v_erL_r + v_erR_r)./2;
    v_eb = zeros(851,2);
    for i = 1:851
        a(i) = cos(vehicleDirection(i));
        b(i) = rad_to_deg * vehicleDirection(i);
        v_eb(i,1) = v_er(1,i) * cos(vehicleDirection(i)); %north
        v_eb(i,2) = v_er(1,i) * sin(vehicleDirection(i)); %east
    end

    position = zeros(851,2);
    for i = 2:851
        position(i,1) = position(i-1,1) + v_eb(i-1,1) * (time(i)-time(i-1)); %north
        position(i,2) = position(i-1,2) + v_eb(i-1,2) * (time(i)-time(i-1)); %east
    end
    
%     figure
%     plot(position(:,2),position(:,1))
%     title('pure INS')
    
    %lat = zeros(851,1);
    %lon = zeros(851,1);
    %lat(1,:) = 51.509253065480780;
    %lon(1,:) = -0.161026458766525;
    %for i = 2:851
        [lat,lon,h] = ned2geodetic(position(:,1)/100000,position(:,2)/100000,zeros(851,1),51.509253065480780,...
            -0.161026458766525,47.114150089444590,referenceEllipsoid);
    %end
    
    figure
    plot(lon,lat);
    title('pure INS navgation');
    xlabel('longitude');
    ylabel('lattidue');
    
    dlmwrite('Pure_INS.csv',[time, lat,lon,v_eb,vehicleDirection],'precision',15);
end