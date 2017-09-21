clear
a = load('GNSS_KF_output.csv');
b = load('Pure_INS.csv');

c = zeros(851,6);
c(:,1) = a(:,1);
c(:,2:6) = 0.65.*a(:,2:6) + 0.35.*b(:,2:6);

figure
plot(c(:,3),c(:,2));
title('INS-GNSS Integrated Navigation');
xlabel('longitude');
ylabel('lattidue');

dlmwrite('INS&GNSS_KF_output.csv',c,'precision',15);