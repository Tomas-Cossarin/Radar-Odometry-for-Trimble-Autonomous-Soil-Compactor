dataset = 1;
IsCalibrated = 'Calibrated';

time_dep_vars = load(strcat("Output_ds"+ dataset +".txt"));
time =      time_dep_vars(:,1);
odo_x =     time_dep_vars(:,2);
test_x_nh = time_dep_vars(:,3);
odo_y =     time_dep_vars(:,4);
test_y_nh = time_dep_vars(:,5);

figure()
plot(time, odo_x, 'r', 'LineWidth', 2)
hold on
plot(time, test_x_nh, 'b', 'LineWidth', 2)
legend('x Position According to Radar Odometry', 'x Position According to GNSS')
ylabel('x Position of the Compactor (m)')
xlabel('Time (s)')
title(strcat('x Position of the Compactor Over Time, Radar Odometry vs GNSS, ', IsCalibrated))

figure()
plot(time, odo_y, 'r', 'LineWidth', 2)
hold on
plot(time, test_y_nh, 'b', 'LineWidth', 2)
legend('y Position According to Radar Odometry', 'y Position According to GNSS')
ylabel('y Position of the Compactor (m)')
xlabel('Time (s)')
title(strcat('y Position of the Compactor Over Time, Radar Odometry vs GNSS, ', IsCalibrated))