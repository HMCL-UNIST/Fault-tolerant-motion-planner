%% Fault Tolerant Motion Planner Results
close all; clear; clc
% Initialize
m = 1.282; g = 9.81;
Ix = 4.856*10^-3; Iy = 4.856*10^-3; Iz = 8.801*10^-3;
b = 0.0066;
kr = 0.02;
kt = 0.5;
l = 0.1;

t_end_mp = 1;
peak_mp = deg2rad(40);
Kp = 10;
Kd = 10;
Kpz = 10;
Kdz = 1;
Ksd = 0.5;
Ks = 10;
Ksdz = 1;
Ksz = 3;
Kb = 4;
Kbz = 1;

sampleRate = 0.001;
numSteps = 100001;
time = sampleRate*[0:(numSteps-1)];
time = time';
tsim = 10;
time1 = 0:sampleRate:tsim;
time1 = time1';

%% Scenario 1-1
time_fault = 1;

des_x = 0; des_y = 0; des_z = 3;
data = [des_x,des_y,des_z].*ones(length(time),1);
des_pos.time = time;
des_pos.signals.values = data;
des_pos.signals.dimensions = 3;
sat = 1;

satInfo = sat*ones(length(time),1);
saturation.time = time;
saturation.signals.values = satInfo;
saturation.signals.dimensions = 1;

phi_transient =30;
phi_info = phi_transient.*ones(length(time),1);
des_phi.time = time;
des_phi.signals.values = phi_info;
des_phi.signals.dimensions = 1;

x = zeros(12,1);
Q.init = [0 0 3 0 0 0  0 0 0 0 0 0]'; %xyz xyzvel pqr phi theta psi
x1 = Q.init(1); x2 = Q.init(2); x3 = Q.init(3);
x4 = Q.init(4); x5 = Q.init(5); x6 = Q.init(6);
x7 = Q.init(7); x8 = Q.init(8); x9 = Q.init(9);
x10 = Q.init(10); x11 = Q.init(11); x12 = Q.init(12);

Q.initInput = [3.1441 3.1441 3.1441 3.1441];
T = m*g;
tp = 0;
tq = 0;
tr = 0;

% Position
x(1) = x4;
x(2) = x5;
x(3) = x6;
% Velocity
x(4) = T/m*(cos(x12)*sin(x11)*cos(x10) + sin(x12)*sin(x10)) - kt*x4;
x(5) = T/m*(sin(x12)*sin(x11)*cos(x10) - cos(x12)*sin(x10)) - kt*x5;
x(6) = -g + T*(1/m)*(cos(x10)*cos(x11)) - kt*x6;
% Body rate
x(7) = ((Ix-Iz)*x8*x9 + tp)/Ix;
x(8) = ((Iz-Ix)*x7*x9 + tq)/Ix;
x(9) = (-kr*x9 + tr)/Iz;
% Angle
x(10) = x7 + x8*sin(x10)*tan(x11) + x9*cos(x10)*tan(x11);
x(11) = x8*cos(x10) - x9*sin(x10);
x(12) = x8*sin(x10)/cos(x11) + x9*cos(x10)*cos(x11);  

Q.initDerivate = x;
saturation_on = 1;

% Pre-failure
[PFTC_state,PFTC_state_dot,PFTC_thrust,pre_err_norm,pre_err_norm_pos,pre_failure_time] = preFailure_func_scen1;

% Update
PFTC_phi = PFTC_state(:,10); PFTC_theta = PFTC_state(:,11); PFTC_psi = PFTC_state(:,12);
PFTC_x = PFTC_state(:,1); PFTC_y = PFTC_state(:,2); PFTC_z = PFTC_state(:,3);

Q.init = PFTC_state(end,:);
Q.initDerivate = PFTC_state_dot(end,:);
Q.initInput = PFTC_thrust(end,:);
stop = length(PFTC_thrust(:,1));

phi_deviation = 10*randn(100,1);

% PD-based FTC
for MC = 1:length(phi_deviation)

    wait = sprintf('Monte Carlo %d/100',MC);
    disp(wait)

    Q.init = [0 0 3 0 0 0  0 0 0 deg2rad(20+phi_deviation(MC)) 0 0]'; %xyz xyzvel pqr phi theta psi
    Q.initInput = [3.1441 3.1441 3.1441 3.1441]/cos(deg2rad(20+phi_deviation(MC)));

    simulation = sim('./Scen1/scen1_PD');
    time_sim = time_sim + time_fault;

    for i = 1:length(x)
        err_norm(MC,i) = norm([phi(i)-deg2rad(0), theta(i)-deg2rad(0)]);
    end
    time_sim = vertcat(time_fault-sampleRate, time_sim);
    F1 = vertcat(PFTC_thrust(end,1),F1);
    F2 = vertcat(PFTC_thrust(end,2),F2);

    simulation = sim('./Scen1/scen1_PDproposed');
    time_sim = time_sim + time_fault;
    for i = 1:length(phi)
        err_norm_proposed(MC,i) = norm([phi(i)-peak_mp, theta(i)-deg2rad(0)]);
    end
    time_sim = vertcat(time_fault-sampleRate, time_sim);
    F1 = vertcat(PFTC_thrust(end,1),F1);
    F2 = vertcat(PFTC_thrust(end,2),F2);
end

err_norm_mean = mean(err_norm);
err_norm_proposed_mean = mean(err_norm_proposed);
err_norm_std = std(err_norm);
err_norm_proposed_std = std(err_norm_proposed);

err_norm_mean = horzcat(err_norm_mean(1),err_norm_mean);
err_norm_proposed_mean = horzcat(err_norm_proposed_mean(1),err_norm_proposed_mean);
err_norm_std = horzcat(err_norm_std(1),err_norm_std);
err_norm_proposed_std = horzcat(err_norm_proposed_std(1),err_norm_proposed_std);

y = rad2deg(err_norm_mean);
x = time_sim';
curve1 = y + rad2deg(err_norm_std);
curve2 = y - rad2deg(err_norm_std);
curve2(curve2 < 0) = 0;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

y_proposed = rad2deg(err_norm_proposed_mean);
curve1_proposed = y_proposed + rad2deg(err_norm_proposed_std);
curve2_proposed = y_proposed - rad2deg(err_norm_proposed_std);
inBetween_proposed = [curve1_proposed, fliplr(curve2_proposed)];

figure(1); hold on; grid on
plot([time_sim(1);time_sim],rad2deg([0 err_norm_mean]),'r--','LineWidth',1)
plot([time_sim(1);time_sim],rad2deg([0 err_norm_proposed_mean]),'b-.','LineWidth',1)
fill(x2, inBetween, 'm', 'FaceAlpha',0.2,'LineStyle','none'); 
fill(x2, inBetween_proposed, 'c', 'FaceAlpha',0.2,'LineStyle','none');
plot(pre_failure_time,rad2deg(pre_err_norm),'k-','LineWidth',1);
legend('Mean (PD)','Mean (PD with proposed motion planner)','STD (PD)','STD (PD with proposed motion planner)')
xlabel('Time [sec]')
ylabel('Attitude error norm [deg]')
ylim([-5 35])
fig1 = figure(1);
fig_size = [500, 300, 600, 350];
set(fig1, 'OuterPosition', fig_size)
xlim([0 6])

% SMC-based FTC
for MC = 1:length(phi_deviation)

    wait = sprintf('Monte Carlo %d/100',MC);
    disp(wait)

    Q.init = [0 0 3 0 0 0  0 0 0 deg2rad(20+phi_deviation(MC)) 0 0]'; %xyz xyzvel pqr phi theta psi
    Q.initInput = [3.1441 3.1441 3.1441 3.1441]/cos(deg2rad(20+phi_deviation(MC)));

    simulation = sim('./Scen1/scen1_SMC');
    time_sim = time_sim + time_fault;

    for i = 1:length(x)
        err_norm(MC,i) = norm([phi(i)-deg2rad(0), theta(i)-deg2rad(0)]);
    end
    time_sim = vertcat(time_fault-sampleRate, time_sim);
    F1 = vertcat(PFTC_thrust(end,1),F1);
    F2 = vertcat(PFTC_thrust(end,2),F2);

    simulation = sim('./Scen1/scen1_SMCproposed');
    time_sim = time_sim + time_fault;
    for i = 1:length(phi)
        err_norm_proposed(MC,i) = norm([phi(i)-peak_mp, theta(i)-deg2rad(0)]);
    end
    time_sim = vertcat(time_fault-sampleRate, time_sim);
    F1 = vertcat(PFTC_thrust(end,1),F1);
    F2 = vertcat(PFTC_thrust(end,2),F2);
end

err_norm_mean = mean(err_norm);
err_norm_proposed_mean = mean(err_norm_proposed);
err_norm_std = std(err_norm);
err_norm_proposed_std = std(err_norm_proposed);

err_norm_mean = horzcat(err_norm_mean(1),err_norm_mean);
err_norm_proposed_mean = horzcat(err_norm_proposed_mean(1),err_norm_proposed_mean);
err_norm_std = horzcat(err_norm_std(1),err_norm_std);
err_norm_proposed_std = horzcat(err_norm_proposed_std(1),err_norm_proposed_std);

y = rad2deg(err_norm_mean);
x = time_sim';
curve1 = y + rad2deg(err_norm_std);
curve2 = y - rad2deg(err_norm_std);
curve2(curve2 < 0) = 0;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

y_proposed = rad2deg(err_norm_proposed_mean);
curve1_proposed = y_proposed + rad2deg(err_norm_proposed_std);
curve2_proposed = y_proposed - rad2deg(err_norm_proposed_std);
inBetween_proposed = [curve1_proposed, fliplr(curve2_proposed)];

figure(2); hold on; grid on
plot([time_sim(1);time_sim],rad2deg([0 err_norm_mean]),'r--','LineWidth',1)
plot([time_sim(1);time_sim],rad2deg([0 err_norm_proposed_mean]),'b-.','LineWidth',1)
fill(x2, inBetween, 'm', 'FaceAlpha',0.2,'LineStyle','none'); 
fill(x2, inBetween_proposed, 'c', 'FaceAlpha',0.2,'LineStyle','none');
plot(pre_failure_time,rad2deg(pre_err_norm),'k-','LineWidth',1);
legend('Mean (SMC)','Mean (SMC with proposed motion planner)','STD (SMC)','STD (SMC with proposed motion planner)')
xlabel('Time [sec]')
ylabel('Attitude error norm [deg]')
ylim([-5 35])
fig2 = figure(2);
fig_size = [500, 300, 600, 350];
set(fig2, 'OuterPosition', fig_size)
xlim([0 6])

% BSC-based FTC
for MC = 1:length(phi_deviation)

    wait = sprintf('Monte Carlo %d/100',MC);
    disp(wait)

    Q.init = [0 0 3 0 0 0  0 0 0 deg2rad(20+phi_deviation(MC)) 0 0]'; %xyz xyzvel pqr phi theta psi
    Q.initInput = [3.1441 3.1441 3.1441 3.1441]/cos(deg2rad(20+phi_deviation(MC)));

    simulation = sim('./Scen1/scen1_BSC');
    time_sim = time_sim + time_fault;

    for i = 1:length(x)
        err_norm(MC,i) = norm([phi(i)-deg2rad(0), theta(i)-deg2rad(0)]);
    end
    time_sim = vertcat(time_fault-sampleRate, time_sim);
    F1 = vertcat(PFTC_thrust(end,1),F1);
    F2 = vertcat(PFTC_thrust(end,2),F2);

    simulation = sim('./Scen1/scen1_BSCproposed');
    time_sim = time_sim + time_fault;
    for i = 1:length(phi)
        err_norm_proposed(MC,i) = norm([phi(i)-peak_mp, theta(i)-deg2rad(0)]);
    end
    time_sim = vertcat(time_fault-sampleRate, time_sim);
    F1 = vertcat(PFTC_thrust(end,1),F1);
    F2 = vertcat(PFTC_thrust(end,2),F2);
end

err_norm_mean = mean(err_norm);
err_norm_proposed_mean = mean(err_norm_proposed);
err_norm_std = std(err_norm);
err_norm_proposed_std = std(err_norm_proposed);

err_norm_mean = horzcat(err_norm_mean(1),err_norm_mean);
err_norm_proposed_mean = horzcat(err_norm_proposed_mean(1),err_norm_proposed_mean);
err_norm_std = horzcat(err_norm_std(1),err_norm_std);
err_norm_proposed_std = horzcat(err_norm_proposed_std(1),err_norm_proposed_std);

y = rad2deg(err_norm_mean);
x = time_sim';
curve1 = y + rad2deg(err_norm_std);
curve2 = y - rad2deg(err_norm_std);
curve2(curve2 < 0) = 0;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

y_proposed = rad2deg(err_norm_proposed_mean);
curve1_proposed = y_proposed + rad2deg(err_norm_proposed_std);
curve2_proposed = y_proposed - rad2deg(err_norm_proposed_std);
inBetween_proposed = [curve1_proposed, fliplr(curve2_proposed)];

figure(3); hold on; grid on
plot([time_sim(1);time_sim],rad2deg([0 err_norm_mean]),'r--','LineWidth',1)
plot([time_sim(1);time_sim],rad2deg([0 err_norm_proposed_mean]),'b-.','LineWidth',1)
fill(x2, inBetween, 'm', 'FaceAlpha',0.2,'LineStyle','none'); 
fill(x2, inBetween_proposed, 'c', 'FaceAlpha',0.2,'LineStyle','none');
plot(pre_failure_time,rad2deg(pre_err_norm),'k-','LineWidth',1);
legend('Mean (BSC)','Mean (BSC with proposed motion planner)','STD (BSC)','STD (BSC with proposed motion planner)')
xlabel('Time [sec]')
ylabel('Attitude error norm [deg]')
ylim([-5 35])
fig3 = figure(3);
fig_size = [500, 300, 600, 350];
set(fig3, 'OuterPosition', fig_size)
xlim([0 6])

%% Scenario 1-2
load('./Data/data_PD_0');
err_norm_PD_max_PD = err_norm_PD_max;
err_z_PD_box_PD = err_z_PD_box;
err_norm_PDproposed_max_PD0 = err_norm_PDproposed_max;
err_z_PDproposed_box_PD0 = err_z_PDproposed_box;

saturation_time_PDall = saturation_time_PD;
saturation_time_PDproposed0 = saturation_time_PDproposed;

load('./Data/data_PD_20');
err_norm_PDproposed_max_PD20 = err_norm_PDproposed_max;
err_z_PDproposed_box_PD20 = err_z_PDproposed_box;

saturation_time_PDproposed20 = saturation_time_PDproposed;

load('./Data/data_PD_40');
err_norm_PDproposed_max_PD40 = err_norm_PDproposed_max;
err_z_PDproposed_box_PD40 = err_z_PDproposed_box;

saturation_time_PDproposed40 = saturation_time_PDproposed;

load('./Data/data_PD_60');
err_norm_PDproposed_max_PD60 = err_norm_PDproposed_max;
err_z_PDproposed_box_PD60 = err_z_PDproposed_box;

saturation_time_PDproposed60 = saturation_time_PDproposed;

load('./Data/data_SMC_0');
err_norm_PD_max_SMC = err_norm_PD_max;
err_z_PD_box_SMC = err_z_PD_box;
err_norm_PDproposed_max_SMC0 = err_norm_PDproposed_max;
err_z_PDproposed_box_SMC0 = err_z_PDproposed_box;

saturation_time_SMCall = saturation_time_PD;
saturation_time_SMCproposed0 = saturation_time_PDproposed;

load('./Data/data_SMC_20');
err_norm_PDproposed_max_SMC20 = err_norm_PDproposed_max;
err_z_PDproposed_box_SMC20 = err_z_PDproposed_box;

saturation_time_SMCproposed20 = saturation_time_PDproposed;

load('./Data/data_SMC_40');
err_norm_PDproposed_max_SMC40 = err_norm_PDproposed_max;
err_z_PDproposed_box_SMC40 = err_z_PDproposed_box;

saturation_time_SMCproposed40 = saturation_time_PDproposed;

load('./Data/data_SMC_60');
err_norm_PDproposed_max_SMC60 = err_norm_PDproposed_max;
err_z_PDproposed_box_SMC60 = err_z_PDproposed_box;

saturation_time_SMCproposed60 = saturation_time_PDproposed;

load('./Data/data_BSC_0');
err_norm_PD_max_FLC = err_norm_PD_max;
err_z_PD_box_FLC = err_z_PD_box;
err_norm_PDproposed_max_FLC0 = err_norm_PDproposed_max;
err_z_PDproposed_box_FLC0 = err_z_PDproposed_box;

saturation_time_FLCall = saturation_time_PD;
saturation_time_FLCproposed0 = saturation_time_PDproposed;

load('./Data/data_BSC_20');
err_norm_PDproposed_max_FLC20 = err_norm_PDproposed_max;
err_z_PDproposed_box_FLC20 = err_z_PDproposed_box;

saturation_time_FLCproposed20 = saturation_time_PDproposed;

load('./Data/data_BSC_40');
err_norm_PDproposed_max_FLC40 = err_norm_PDproposed_max;
err_z_PDproposed_box_FLC40 = err_z_PDproposed_box;

saturation_time_FLCproposed40 = saturation_time_PDproposed;

load('./Data/data_BSC_60');
err_norm_PDproposed_max_FLC60 = err_norm_PDproposed_max;
err_z_PDproposed_box_FLC60 = err_z_PDproposed_box;

saturation_time_FLCproposed60 = saturation_time_PDproposed;

cmap = hsv;

% Box plot
err_norm_PDproposed_max_PDall = min(min(min(err_norm_PDproposed_max_PD0,err_norm_PDproposed_max_PD20),err_norm_PDproposed_max_PD40),err_norm_PDproposed_max_PD60);
err_norm_PDproposed_max_SMCall = min(min(min(err_norm_PDproposed_max_SMC0,err_norm_PDproposed_max_SMC20),err_norm_PDproposed_max_SMC40),err_norm_PDproposed_max_SMC60);
err_norm_PDproposed_max_FLCall = min(min(min(err_norm_PDproposed_max_FLC0,err_norm_PDproposed_max_FLC20),err_norm_PDproposed_max_FLC40),err_norm_PDproposed_max_FLC60);
min1_PD = err_z_PDproposed_box_PD0 .* (abs(err_z_PDproposed_box_PD0) <= abs(err_z_PDproposed_box_PD20)) + err_z_PDproposed_box_PD20 .* (abs(err_z_PDproposed_box_PD20) < abs(err_z_PDproposed_box_PD0));
min2_PD = min1_PD .* (abs(min1_PD) <= abs(err_z_PDproposed_box_PD40)) + err_z_PDproposed_box_PD40 .* (abs(err_z_PDproposed_box_PD40) < abs(min1_PD));
err_z_PDproposed_box_PDall = min2_PD .* (abs(min2_PD) <= abs(err_z_PDproposed_box_PD60)) + err_z_PDproposed_box_PD60 .* (abs(err_z_PDproposed_box_PD60) < abs(min2_PD));
min1_SMC = err_z_PDproposed_box_SMC0 .* (abs(err_z_PDproposed_box_SMC0) <= abs(err_z_PDproposed_box_SMC20)) + err_z_PDproposed_box_SMC20 .* (abs(err_z_PDproposed_box_SMC20) < abs(err_z_PDproposed_box_SMC0));
min2_SMC = min1_SMC .* (abs(min1_SMC) <= abs(err_z_PDproposed_box_SMC40)) + err_z_PDproposed_box_SMC40 .* (abs(err_z_PDproposed_box_SMC40) < abs(min1_SMC));
err_z_PDproposed_box_SMCall = min2_SMC .* (abs(min2_SMC) <= abs(err_z_PDproposed_box_SMC60)) + err_z_PDproposed_box_SMC60 .* (abs(err_z_PDproposed_box_SMC60) < abs(min2_SMC));
min1_FLC = err_z_PDproposed_box_FLC0 .* (abs(err_z_PDproposed_box_FLC0) <= abs(err_z_PDproposed_box_FLC20)) + err_z_PDproposed_box_FLC20 .* (abs(err_z_PDproposed_box_FLC20) < abs(err_z_PDproposed_box_FLC0));
min2_FLC = min1_FLC .* (abs(min1_FLC) <= abs(err_z_PDproposed_box_FLC40)) + err_z_PDproposed_box_FLC40 .* (abs(err_z_PDproposed_box_FLC40) < abs(min1_FLC));
err_z_PDproposed_box_FLCall = min2_FLC .* (abs(min2_FLC) <= abs(err_z_PDproposed_box_FLC60)) + err_z_PDproposed_box_FLC60 .* (abs(err_z_PDproposed_box_FLC60) < abs(min2_FLC));

i = 6;
x = 1:i;
colors = [[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];[0.8500 0.3250 0.0980];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0 0.4470 0.7410]];

figure, boxplot([rad2deg(err_norm_PD_max_PD(:)) rad2deg(err_norm_PDproposed_max_PDall(:)) rad2deg(err_norm_PD_max_SMC(:)) rad2deg(err_norm_PDproposed_max_SMCall(:)) rad2deg(err_norm_PD_max_FLC(:)) rad2deg(err_norm_PDproposed_max_FLCall(:))],'whisker',inf,'Labels',{'PD','Proposed PD','SMC','Proposed SMC','BSC','Proposed BSC'}), ylabel('Maximum attitude error norm [deg]')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

figure, boxplot([err_z_PD_box_PD(:) err_z_PDproposed_box_PDall(:) err_z_PD_box_SMC(:) err_z_PDproposed_box_SMCall(:) err_z_PD_box_FLC(:) err_z_PDproposed_box_FLCall(:)],'whisker',inf,'Labels',{'PD','Proposed PD','SMC','Proposed SMC','BSC','Proposed BSC'}), ylabel('Maximum altitude error [m]')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

% Bar graph
saturation_time_PDproposedall = min(min(min(saturation_time_PDproposed0,saturation_time_PDproposed20),saturation_time_PDproposed40),saturation_time_PDproposed60);
saturation_time_SMCproposedall = min(min(min(saturation_time_SMCproposed0,saturation_time_SMCproposed20),saturation_time_SMCproposed40),saturation_time_SMCproposed60);
saturation_time_FLCproposedall = min(min(min(saturation_time_FLCproposed0,saturation_time_FLCproposed20),saturation_time_FLCproposed40),saturation_time_FLCproposed60);

stackData = [];
stackData(1,1,1) = sum(saturation_time_PDall(:))/(size(saturation_time_PDall,1)*size(saturation_time_PDall,2));
stackData(1,2,1) = sum(saturation_time_PDproposedall(:))/(size(saturation_time_PDproposedall,1)*size(saturation_time_PDproposedall,2));
stackData(2,1,1) = sum(saturation_time_SMCall(:))/(size(saturation_time_SMCall,1)*size(saturation_time_SMCall,2));
stackData(2,2,1) = sum(saturation_time_SMCproposedall(:))/(size(saturation_time_SMCproposedall,1)*size(saturation_time_SMCproposedall,2));
stackData(3,1,1) = sum(saturation_time_FLCall(:))/(size(saturation_time_FLCall,1)*size(saturation_time_FLCall,2));
stackData(3,2,1) = sum(saturation_time_FLCproposedall(:))/(size(saturation_time_FLCproposedall,1)*size(saturation_time_FLCproposedall,2));

plotBarStackGroups(stackData/15*100,{'PD';'SMC';'BSC'});
ylabel('Average time ratio of FTC input being saturated [%]'), legend('Conventional FTC','FTC with proposed motion planner')

% heat map
load('./Data/data_PD_0');
err_norm_PDproposed_max0 = err_norm_PDproposed_max;
err_z_PDproposed_box0 = err_z_PDproposed_box;
safe_PDproposed_max0 = safe_PDproposed_max;
safe_z_PDproposed_max0 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max0 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max0 = safe_saturation_binary_PDproposed_max;

load('./Data/data_PD_20');
err_norm_PDproposed_max20 = err_norm_PDproposed_max;
err_z_PDproposed_box20 = err_z_PDproposed_box;
safe_PDproposed_max20 = safe_PDproposed_max;
safe_z_PDproposed_max20 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max20 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max20 = safe_saturation_binary_PDproposed_max;

load('./Data/data_PD_40');
err_norm_PDproposed_max40 = err_norm_PDproposed_max;
err_z_PDproposed_box40 = err_z_PDproposed_box;
safe_PDproposed_max40 = safe_PDproposed_max;
safe_z_PDproposed_max40 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max40 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max40 = safe_saturation_binary_PDproposed_max;

load('./Data/data_PD_60');
err_norm_PDproposed_max60 = err_norm_PDproposed_max;
err_z_PDproposed_box60 = err_z_PDproposed_box;
safe_PDproposed_max60 = safe_PDproposed_max;
safe_z_PDproposed_max60 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max60 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max60 = safe_saturation_binary_PDproposed_max;

cmap = hsv;

for i = 1:size(safe_PDproposed_max0,1)
    for j = 1:size(safe_PDproposed_max0,2)
        if safe_PDproposed_max0(i,j) == 200
            safe_PDproposed_max0(i,j) = 0;
            safe_z_PDproposed_max0(i,j) = 0;
            safe_saturation_PDproposed_max0(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max20,1)
    for j = 1:size(safe_PDproposed_max20,2)
        if safe_PDproposed_max20(i,j) == 200
            safe_PDproposed_max20(i,j) = 0;
            safe_z_PDproposed_max20(i,j) = 0;
            safe_saturation_PDproposed_max20(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max40,1)
    for j = 1:size(safe_PDproposed_max40,2)
        if safe_PDproposed_max40(i,j) == 200
            safe_PDproposed_max40(i,j) = 0;
            safe_z_PDproposed_max40(i,j) = 0;
            safe_saturation_PDproposed_max40(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max60,1)
    for j = 1:size(safe_PDproposed_max60,2)
        if safe_PDproposed_max60(i,j) == 200
            safe_PDproposed_max60(i,j) = 0;
            safe_z_PDproposed_max60(i,j) = 0;
            safe_saturation_PDproposed_max60(i,j) = 0;
        end
    end
end

err_norm_PDproposed_max_all = min(min(min(err_norm_PDproposed_max0,err_norm_PDproposed_max20),err_norm_PDproposed_max40),err_norm_PDproposed_max60);
err_z_PDproposed_box_all = min(min(min(err_z_PDproposed_box0,err_z_PDproposed_box20),err_z_PDproposed_box40),err_z_PDproposed_box60);
safe_PDproposed_max_all = max(max(max(safe_PDproposed_max0,safe_PDproposed_max20),safe_PDproposed_max40),safe_PDproposed_max60);
safe_z_PDproposed_max_all = max(max(max(safe_z_PDproposed_max0,safe_z_PDproposed_max20),safe_z_PDproposed_max40),safe_z_PDproposed_max60);
safe_saturation_PDproposed_max_all = max(max(max(safe_saturation_PDproposed_max0,safe_saturation_PDproposed_max20),safe_saturation_PDproposed_max40),safe_saturation_PDproposed_max60);
safe_saturation_binary_PDproposed_max_all = max(max(max(safe_saturation_binary_PDproposed_max0,safe_saturation_binary_PDproposed_max20),safe_saturation_binary_PDproposed_max40),safe_saturation_binary_PDproposed_max60);

A = safe_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('PD attitude error norm')

A = safe_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('PD proposed attitude error norm')

A = safe_z_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('PD altitude error')

A = safe_z_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('PD proposed altitude error')

A = safe_saturation_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('PD saturation')

A = safe_saturation_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('PD proposed saturation')

load('./Data/data_SMC_0');
err_norm_PDproposed_max0 = err_norm_PDproposed_max;
err_z_PDproposed_box0 = err_z_PDproposed_box;
safe_PDproposed_max0 = safe_PDproposed_max;
safe_z_PDproposed_max0 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max0 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max0 = safe_saturation_binary_PDproposed_max;

load('./Data/data_SMC_20');
err_norm_PDproposed_max20 = err_norm_PDproposed_max;
err_z_PDproposed_box20 = err_z_PDproposed_box;
safe_PDproposed_max20 = safe_PDproposed_max;
safe_z_PDproposed_max20 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max20 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max20 = safe_saturation_binary_PDproposed_max;

load('./Data/data_SMC_40');
err_norm_PDproposed_max40 = err_norm_PDproposed_max;
err_z_PDproposed_box40 = err_z_PDproposed_box;
safe_PDproposed_max40 = safe_PDproposed_max;
safe_z_PDproposed_max40 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max40 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max40 = safe_saturation_binary_PDproposed_max;

load('./Data/data_SMC_60');
err_norm_PDproposed_max60 = err_norm_PDproposed_max;
err_z_PDproposed_box60 = err_z_PDproposed_box;
safe_PDproposed_max60 = safe_PDproposed_max;
safe_z_PDproposed_max60 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max60 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max60 = safe_saturation_binary_PDproposed_max;

cmap = hsv;

for i = 1:size(safe_PDproposed_max0,1)
    for j = 1:size(safe_PDproposed_max0,2)
        if safe_PDproposed_max0(i,j) == 200
            safe_PDproposed_max0(i,j) = 0;
            safe_z_PDproposed_max0(i,j) = 0;
            safe_saturation_PDproposed_max0(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max20,1)
    for j = 1:size(safe_PDproposed_max20,2)
        if safe_PDproposed_max20(i,j) == 200
            safe_PDproposed_max20(i,j) = 0;
            safe_z_PDproposed_max20(i,j) = 0;
            safe_saturation_PDproposed_max20(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max40,1)
    for j = 1:size(safe_PDproposed_max40,2)
        if safe_PDproposed_max40(i,j) == 200
            safe_PDproposed_max40(i,j) = 0;
            safe_z_PDproposed_max40(i,j) = 0;
            safe_saturation_PDproposed_max40(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max60,1)
    for j = 1:size(safe_PDproposed_max60,2)
        if safe_PDproposed_max60(i,j) == 200
            safe_PDproposed_max60(i,j) = 0;
            safe_z_PDproposed_max60(i,j) = 0;
            safe_saturation_PDproposed_max60(i,j) = 0;
        end
    end
end

err_norm_PDproposed_max_all = min(min(min(err_norm_PDproposed_max0,err_norm_PDproposed_max20),err_norm_PDproposed_max40),err_norm_PDproposed_max60);
err_z_PDproposed_box_all = min(min(min(err_z_PDproposed_box0,err_z_PDproposed_box20),err_z_PDproposed_box40),err_z_PDproposed_box60);
safe_PDproposed_max_all = max(max(max(safe_PDproposed_max0,safe_PDproposed_max20),safe_PDproposed_max40),safe_PDproposed_max60);
safe_z_PDproposed_max_all = max(max(max(safe_z_PDproposed_max0,safe_z_PDproposed_max20),safe_z_PDproposed_max40),safe_z_PDproposed_max60);
safe_saturation_PDproposed_max_all = max(max(max(safe_saturation_PDproposed_max0,safe_saturation_PDproposed_max20),safe_saturation_PDproposed_max40),safe_saturation_PDproposed_max60);
safe_saturation_binary_PDproposed_max_all = max(max(max(safe_saturation_binary_PDproposed_max0,safe_saturation_binary_PDproposed_max20),safe_saturation_binary_PDproposed_max40),safe_saturation_binary_PDproposed_max60);

A = safe_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('SMC attitude error norm')

A = safe_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('SMC proposed attitude error norm')

A = safe_z_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('SMC altitude error')

A = safe_z_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('SMC proposed altitude error')

A = safe_saturation_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('SMC saturation')

A = safe_saturation_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('SMC proposed saturation')

load('./Data/data_BSC_0');
err_norm_PDproposed_max0 = err_norm_PDproposed_max;
err_z_PDproposed_box0 = err_z_PDproposed_box;
safe_PDproposed_max0 = safe_PDproposed_max;
safe_z_PDproposed_max0 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max0 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max0 = safe_saturation_binary_PDproposed_max;

load('./Data/data_BSC_20');
err_norm_PDproposed_max20 = err_norm_PDproposed_max;
err_z_PDproposed_box20 = err_z_PDproposed_box;
safe_PDproposed_max20 = safe_PDproposed_max;
safe_z_PDproposed_max20 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max20 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max20 = safe_saturation_binary_PDproposed_max;

load('./Data/data_BSC_40');
err_norm_PDproposed_max40 = err_norm_PDproposed_max;
err_z_PDproposed_box40 = err_z_PDproposed_box;
safe_PDproposed_max40 = safe_PDproposed_max;
safe_z_PDproposed_max40 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max40 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max40 = safe_saturation_binary_PDproposed_max;

load('./Data/data_BSC_60');
err_norm_PDproposed_max60 = err_norm_PDproposed_max;
err_z_PDproposed_box60 = err_z_PDproposed_box;
safe_PDproposed_max60 = safe_PDproposed_max;
safe_z_PDproposed_max60 = safe_z_PDproposed_max;
safe_saturation_PDproposed_max60 = safe_saturation_PDproposed_max;
safe_saturation_binary_PDproposed_max60 = safe_saturation_binary_PDproposed_max;

cmap = hsv;

for i = 1:size(safe_PDproposed_max0,1)
    for j = 1:size(safe_PDproposed_max0,2)
        if safe_PDproposed_max0(i,j) == 200
            safe_PDproposed_max0(i,j) = 0;
            safe_z_PDproposed_max0(i,j) = 0;
            safe_saturation_PDproposed_max0(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max20,1)
    for j = 1:size(safe_PDproposed_max20,2)
        if safe_PDproposed_max20(i,j) == 200
            safe_PDproposed_max20(i,j) = 0;
            safe_z_PDproposed_max20(i,j) = 0;
            safe_saturation_PDproposed_max20(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max40,1)
    for j = 1:size(safe_PDproposed_max40,2)
        if safe_PDproposed_max40(i,j) == 200
            safe_PDproposed_max40(i,j) = 0;
            safe_z_PDproposed_max40(i,j) = 0;
            safe_saturation_PDproposed_max40(i,j) = 0;
        end
    end
end
for i = 1:size(safe_PDproposed_max60,1)
    for j = 1:size(safe_PDproposed_max60,2)
        if safe_PDproposed_max60(i,j) == 200
            safe_PDproposed_max60(i,j) = 0;
            safe_z_PDproposed_max60(i,j) = 0;
            safe_saturation_PDproposed_max60(i,j) = 0;
        end
    end
end

err_norm_PDproposed_max_all = min(min(min(err_norm_PDproposed_max0,err_norm_PDproposed_max20),err_norm_PDproposed_max40),err_norm_PDproposed_max60);
err_z_PDproposed_box_all = min(min(min(err_z_PDproposed_box0,err_z_PDproposed_box20),err_z_PDproposed_box40),err_z_PDproposed_box60);
safe_PDproposed_max_all = max(max(max(safe_PDproposed_max0,safe_PDproposed_max20),safe_PDproposed_max40),safe_PDproposed_max60);
safe_z_PDproposed_max_all = max(max(max(safe_z_PDproposed_max0,safe_z_PDproposed_max20),safe_z_PDproposed_max40),safe_z_PDproposed_max60);
safe_saturation_PDproposed_max_all = max(max(max(safe_saturation_PDproposed_max0,safe_saturation_PDproposed_max20),safe_saturation_PDproposed_max40),safe_saturation_PDproposed_max60);
safe_saturation_binary_PDproposed_max_all = max(max(max(safe_saturation_binary_PDproposed_max0,safe_saturation_binary_PDproposed_max20),safe_saturation_binary_PDproposed_max40),safe_saturation_binary_PDproposed_max60);

A = safe_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('BSC attitude error norm')

A = safe_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('BSC proposed attitude error norm')

A = safe_z_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('BSC altitude error')

A = safe_z_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('BSC proposed altitude error')

A = safe_saturation_PD_max;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
xlabel('Initial Pitch [deg]'), ylabel('Initial Roll [deg]')
axis square
title('BSC saturation')

A = safe_saturation_PDproposed_max_all;
figure
image(A)
colormap(cmap)
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','normal')
axis square
title('BSC proposed saturation')


%% Scenario 2
clear all
m = 1.282; g = 9.81;
Ix = 4.856*10^-3; Iy = 4.856*10^-3; Iz = 8.801*10^-3;
b = 0.0066;
kr = 0.02;
kt = 0.5;
l = 0.1;

t_end_mp = 1;
peak_mp = deg2rad(40);
Kpxy = 0.1;
Kp = 10;
Kd = 10;
Kpz = 10;
Kdz = 1;
Ksd = 0.5;
Ks = 10;
Ksdz = 1;
Ksz = 3;
Kb = 4;
Kbz = 1;

sampleRate = 0.001;
numSteps = 100001;
time = sampleRate*[0:(numSteps-1)];
time = time';
tsim = 10;
time1 = 0:sampleRate:tsim;
time1 = time1';

time_fault = 5;

% desired position
init_x = 0; init_y = 0; init_z = 3;
des_x = 0; des_y = 0; des_z = 3;
data = [des_x,des_y,des_z].*ones(length(time),1);
des_pos.time = time;
des_pos.signals.values = data;
des_pos.signals.dimensions = 3;
sat = 1; %Set saturation  == 1, No saturation == 0;

satInfo = sat*ones(length(time),1);
saturation.time = time;
saturation.signals.values = satInfo;
saturation.signals.dimensions = 1;

phi_transient =30;
phi_info = phi_transient.*ones(length(time),1);
des_phi.time = time;
des_phi.signals.values = phi_info;
des_phi.signals.dimensions = 1;


x = zeros(12,1);
Q.init = [0 0 3 0 0 0  0 0 0 0 0 0]'; %xyz xyzvel pqr phi theta psi
x1 = Q.init(1); x2 = Q.init(2); x3 = Q.init(3);
x4 = Q.init(4); x5 = Q.init(5); x6 = Q.init(6);
x7 = Q.init(7); x8 = Q.init(8); x9 = Q.init(9);
x10 = Q.init(10); x11 = Q.init(11); x12 = Q.init(12);

Q.initInput = [3.1441 3.1441 3.1441 3.1441];
T = m*g;
tp = 0;
tq = 0;
tr = 0;

% Position
x(1) = x4;
x(2) = x5;
x(3) = x6;
% Velocity
x(4) = T/m*(cos(x12)*sin(x11)*cos(x10) + sin(x12)*sin(x10)) - kt*x4;
x(5) = T/m*(sin(x12)*sin(x11)*cos(x10) - cos(x12)*sin(x10)) - kt*x5;
x(6) = -g + T*(1/m)*(cos(x10)*cos(x11)) - kt*x6;
% Body rate
x(7) = ((Ix-Iz)*x8*x9 + tp)/Ix;
x(8) = ((Iz-Ix)*x7*x9 + tq)/Ix;
x(9) = (-kr*x9 + tr)/Iz;
% Angle
x(10) = x7 + x8*sin(x10)*tan(x11) + x9*cos(x10)*tan(x11);
x(11) = x8*cos(x10) - x9*sin(x10);
x(12) = x8*sin(x10)/cos(x11) + x9*cos(x10)*cos(x11);  

Q.initDerivate = x;
saturation_on = 1;

% Pre-failure
[PFTC_state,PFTC_state_dot,PFTC_thrust,pre_err_norm,pre_err_norm_pos,pre_failure_time,...
 ax25,ax25_enlarge,ax26,ax26_enlarge,ax27,ax27_enlarge] = preFailure_func_scen2;

% Update
PFTC_phi = PFTC_state(:,10); PFTC_theta = PFTC_state(:,11); PFTC_psi = PFTC_state(:,12);
PFTC_x = PFTC_state(:,1); PFTC_y = PFTC_state(:,2); PFTC_z = PFTC_state(:,3);

Q.init = PFTC_state(end,:);
Q.initDerivate = PFTC_state_dot(end,:);
Q.initInput = PFTC_thrust(end,:);
stop = length(PFTC_thrust(:,1));

tsim = 80;
time1 = 0:sampleRate:tsim;
time1 = time1';
des_x = 25; des_y = -10; des_z = 3;
data = [des_x,des_y,des_z].*ones(length(time),1);
des_pos.time = time;
des_pos.signals.values = data;
des_pos.signals.dimensions = 3;

% PD-based FTC
simulation = sim('./Scen2/scen2_PD.slx');
clearvars err_norm
time_sim = time_sim + time_fault;
for i = 1:length(x)
    err_norm(i) = norm([x(i)-des_x, y(i)-des_y, z(i)-des_z]);
end
time_sim = vertcat(time_fault-sampleRate, time_sim);
F1 = vertcat(PFTC_thrust(end,1),F1);
F2 = vertcat(PFTC_thrust(end,2),F2);
err_norm = horzcat(pre_err_norm_pos(end),err_norm);
x = vertcat(PFTC_x(end),x);
y = vertcat(PFTC_y(end),y);
z = vertcat(PFTC_z(end),z);

figure(25); hold on; grid on;
axes(ax25)
plot(time_sim,F1,'r--','LineWidth',1)
plot(time_sim, F2, 'r:','LineWidth',1)
axes(ax25_enlarge)
plot(time_sim,F2,'r:','LineWidth',1.2)

figure(28); hold on
plot3(x,y,z,'b--','LineWidth',1)

x_fall_PD = x(end);
y_fall_PD = y(end);
z_fall_PD = z(end);

alpha = 0.02;
phiPE = peak_mp;
phiPE2 = phiPE - alpha;
rPE2 = (-kr/b+sqrt((kr/b)^2+8*m*g/(l*cos(phiPE2))*(Iz-Iy)*tan(phiPE2)))/(4*(Iz-Iy)*tan(phiPE2)/l);
psidotPE2 = rPE2*(tan(phiPE2)*sin(phiPE2)+cos(phiPE2)); 
TPE = 2*pi/psidotPE2;

% PD-based Proposed method
simulation = sim('./Scen2/scen2_PDproposed');
clearvars err_norm
time_sim = time_sim + time_fault;
for i = 1:length(phi)
    err_norm(i) = norm([x(i)-des_x, y(i)-des_y, z(i)-des_z]);
end
time_sim = vertcat(time_fault-sampleRate, time_sim);
F1 = vertcat(PFTC_thrust(end,1),F1);
F2 = vertcat(PFTC_thrust(end,2),F2);
err_norm = horzcat(pre_err_norm_pos(end),err_norm);
x = vertcat(PFTC_x(end),x);
y = vertcat(PFTC_y(end),y);
z = vertcat(PFTC_z(end),z);

figure(25);
axes(ax25)
plot(time_sim,F1,'b-','LineWidth',1)
plot(time_sim, F2, 'b-.','LineWidth',1)
ylim([-1 14])
fig25 = figure(25);
fig_size = [500, 300, 600, 360];
set(fig25, 'OuterPosition', fig_size)
legend('F1 at Pre-failure','F2 at Pre-failure','F1 of PD','F2 of PD','F1 of PD with Proposed Motion Planner','F2 of PD with Proposed Motion Planner','Location','NorthEast')
ylabel('Rotor thrusts [N]')
xlabel('Time [sec]')
xlim([0 60])

axes(ax25_enlarge)
plot(time_sim,F2,'b-.','LineWidth',1.1)
ylim([-0.1 0.5])
xlim([5.001 6.5])
axes(ax25)
grid on
axes(ax25_enlarge)

figure(28)
plot3(x,y,z,'b','LineWidth',1)

% SM-based FTC
simulation = sim('./Scen2/scen2_SMC');
clearvars err_norm
time_sim = time_sim + time_fault;
for i = 1:length(x)
    err_norm(i) = norm([x(i)-des_x, y(i)-des_y, z(i)-des_z]);
end
time_sim = vertcat(time_fault-sampleRate, time_sim);
F1 = vertcat(PFTC_thrust(end,1),F1);
F2 = vertcat(PFTC_thrust(end,2),F2);
err_norm = horzcat(pre_err_norm_pos(end),err_norm);
x = vertcat(PFTC_x(end),x);
y = vertcat(PFTC_y(end),y);
z = vertcat(PFTC_z(end),z);

figure(26); hold on; grid on;
axes(ax26)
plot(time_sim,F1,'r--','LineWidth',1)
plot(time_sim, F2, 'r:','LineWidth',1)
axes(ax26_enlarge)
plot(time_sim,F2,'r:','LineWidth',1.2)

figure(28)
plot3(x,y,z,'r--','LineWidth',1)

x_fall_SMC = x(end);
y_fall_SMC = y(end);
z_fall_SMC = z(end);

% SM-based Proposed method
simulation = sim('./Scen2/scen2_SMCproposed');
clearvars err_norm
time_sim = time_sim + time_fault;
for i = 1:length(phi)
    err_norm(i) = norm([x(i)-des_x, y(i)-des_y, z(i)-des_z]);
end
time_sim = vertcat(time_fault-sampleRate, time_sim);
F1 = vertcat(PFTC_thrust(end,1),F1);
F2 = vertcat(PFTC_thrust(end,2),F2);
err_norm = horzcat(pre_err_norm_pos(end),err_norm);
x = vertcat(PFTC_x(end),x);
y = vertcat(PFTC_y(end),y);
z = vertcat(PFTC_z(end),z);

figure(26);
axes(ax26)
plot(time_sim,F1,'b-','LineWidth',1)
plot(time_sim, F2, 'b-.','LineWidth',1)
ylim([-1 14])
fig26 = figure(26);
fig_size = [500, 300, 600, 370];
set(fig26, 'OuterPosition', fig_size)
legend('F1 at Pre-failure','F2 at Pre-failure','F1 of SMC','F2 of SMC','F1 of SMC with Proposed Motion Planner','F2 of SMC with Proposed Motion Planner','Location','NorthEast')
ylabel('Rotor thrusts [N]')
xlabel('Time [sec]')
xlim([0 60])

axes(ax26_enlarge)
plot(time_sim,F2,'b-.','LineWidth',1.1)
ylim([-0.1 0.5])
xlim([5.001 6.5])
axes(ax26)
grid on
axes(ax26_enlarge)

figure(28)
plot3(x,y,z,'r','LineWidth',1)

simulation = sim('./Scen2/scen2_BSC');
clearvars err_norm
time_sim = time_sim + time_fault;
for i = 1:length(x)
    err_norm(i) = norm([x(i)-des_x, y(i)-des_y, z(i)-des_z]);
end
time_sim = vertcat(time_fault-sampleRate, time_sim);
F1 = vertcat(PFTC_thrust(end,1),F1);
F2 = vertcat(PFTC_thrust(end,2),F2);
err_norm = horzcat(pre_err_norm_pos(end),err_norm);
x = vertcat(PFTC_x(end),x);
y = vertcat(PFTC_y(end),y);
z = vertcat(PFTC_z(end),z);

figure(27); hold on; grid on;
axes(ax27)
plot(time_sim,F1,'r--','LineWidth',1)
plot(time_sim, F2, 'r:','LineWidth',1)
axes(ax27_enlarge)
plot(time_sim,F2,'r:','LineWidth',1.2)

figure(28)
plot3(x,y,z,'g--','LineWidth',1)

x_fall_BS = x(end);
y_fall_BS = y(end);
z_fall_BS = z(end);

des_position = [des_x des_y];

simulation = sim('./Scen2/scen2_BSCproposed');
clearvars err_norm
time_sim = time_sim + time_fault;
for i = 1:length(phi)
    err_norm(i) = norm([x(i)-des_x, y(i)-des_y, z(i)-des_z]);
end
time_sim = vertcat(time_fault-sampleRate, time_sim);
F1 = vertcat(PFTC_thrust(end,1),F1);
F2 = vertcat(PFTC_thrust(end,2),F2);
err_norm = horzcat(pre_err_norm_pos(end),err_norm);
x = vertcat(PFTC_x(end),x);
y = vertcat(PFTC_y(end),y);
z = vertcat(PFTC_z(end),z);

figure(27);
axes(ax27)
plot(time_sim,F1,'b-','LineWidth',1)
plot(time_sim, F2, 'b-.','LineWidth',1)
fig27 = figure(27);
fig_size = [500, 300, 600, 360];
set(fig27, 'OuterPosition', fig_size)
legend('F1 at Pre-failure','F2 at Pre-failure','F1 of BSC','F2 of BSC','F1 of BSC with Proposed Motion Planner','F2 of BSC with Proposed Motion Planner','Location','NorthEast')
ylabel('Rotor thrusts [N]')
xlabel('Time [sec]')
xlim([0 60])
ylim([-1 14])

axes(ax27_enlarge)
plot(time_sim,F2,'b-.','LineWidth',1.1)
ylim([-0.1 0.5])
xlim([5.001 6.5])
axes(ax27)
grid on
axes(ax27_enlarge)

figure(28)
plot3(x,y,z,'g','LineWidth',1)

r = 1; n = 100;   
[X,Y,Z] = cylinder(r,n);
h = 3;
X = X+des_x;
Y = Y+des_y;
Z = Z*h;

fill3(X(2,:),Y(2,:),Z(2,:),[200/255 1 1],'EdgeColor','k')

plot3(x_fall_PD,y_fall_PD,z_fall_PD,'bx','MarkerSize',12)
plot3(x_fall_SMC,y_fall_SMC,z_fall_SMC,'rx','MarkerSize',12)
plot3(x_fall_BS,y_fall_BS,z_fall_BS,'gx','MarkerSize',12)

xlim([-2 30]), ylim([-15 5])
xlabel('X [m]'), ylabel('Y [m]'), zlabel('Z [m]')
legend('Initial Position','Pre-failure','PD','PD with Proposed Motion Planner','SMC','SMC with Proposed Motion Planner','BSC','BSC with Proposed Motion Planner','Safe zone','location','southwest')

%% Scenario 3
clear all 
close all

m = 1.282; g = 9.81;
Ix = 4.856*10^-3; Iy = 4.856*10^-3; Iz = 8.801*10^-3;
b = 0.0066;
kr = 0.02;
kt = 0.5;
l = 0.1;

t_end_mp = 1;
peak_mp = deg2rad(40);

Kp = 10;
Kd = 10;
Kpz = 10;
Kdz = 1;
Ksd = 0.5;
Ks = 10;
Ksdz = 1;
Ksz = 3;
Kb = 4;
Kbz = 1;
Ab = 22;
Abz = Ab;
Gb = 100;
Gbz = Gb;

sampleRate = 0.01;
numSteps = 100001;
total_time = sampleRate*[0:(numSteps-1)];
total_time = total_time';
tsim = 10;
time1 = 0:sampleRate:tsim;
time1 = time1';

% desired position
init_x = 0; init_y = 0; init_z = 3;
des_x = 0; des_y = 0; des_z = 3;
data = [des_x,des_y,des_z].*ones(length(total_time),1);
des_pos.time = total_time;
des_pos.signals.values = data;
des_pos.signals.dimensions = 3;
sat = 1;

satInfo = sat*ones(length(total_time),1);
saturation.time = total_time;
saturation.signals.values = satInfo;
saturation.signals.dimensions = 1;

phi_transient =30;
phi_info = phi_transient.*ones(length(total_time),1);
des_phi.time = total_time;
des_phi.signals.values = phi_info;
des_phi.signals.dimensions = 1;

x = zeros(12,1);
Q.init = [0 0 3 0 0 0  0 0 0 0 0 0]'; %xyz xyzvel pqr phi theta psi
x1 = Q.init(1); x2 = Q.init(2); x3 = Q.init(3);
x4 = Q.init(4); x5 = Q.init(5); x6 = Q.init(6);
x7 = Q.init(7); x8 = Q.init(8); x9 = Q.init(9);
x10 = Q.init(10); x11 = Q.init(11); x12 = Q.init(12);

Q.initInput = [3.1441 3.1441 3.1441 3.1441];
T = m*g;
tp = 0;
tq = 0;
tr = 0;

% Position
x(1) = x4;
x(2) = x5;
x(3) = x6;
% Velocity
x(4) = T/m*(cos(x12)*sin(x11)*cos(x10) + sin(x12)*sin(x10)) - kt*x4;
x(5) = T/m*(sin(x12)*sin(x11)*cos(x10) - cos(x12)*sin(x10)) - kt*x5;
x(6) = -g + T*(1/m)*(cos(x10)*cos(x11)) - kt*x6;
% Body rate
x(7) = ((Ix-Iz)*x8*x9 + tp)/Ix;
x(8) = ((Iz-Ix)*x7*x9 + tq)/Ix;
x(9) = (-kr*x9 + tr)/Iz;
% Angle
x(10) = x7 + x8*sin(x10)*tan(x11) + x9*cos(x10)*tan(x11);
x(11) = x8*cos(x10) - x9*sin(x10);
x(12) = x8*sin(x10)/cos(x11) + x9*cos(x10)*cos(x11);

Q.initDerivate = x;
saturation_on = 1;

alpha = 0.02;

tsim = 9000;
time1 = 0:sampleRate:tsim;
time1 = time1';
load('trajectory')


index2 = 1;
mu = 0;
sigma = 0;
noise = sigma*randn(1000000,3);


simulation = sim('./Scen3/scen3_PD.slx');
runtime2 = runtime;
filename = sprintf('./Data/reference_%02d.mat', t);
save(filename);

figure(29); hold on
hold on; grid on;
plot(x,y, 'b-')


tsim = 240;
time1 = 0:sampleRate:tsim;
time1 = time1';
load('trajectory')

index2 = 1;
noise = sigma*randn(80021,3);
simulation = sim('./Scen3/scen3_PDproposed');
filename = sprintf('./Data/proposed_%02d.mat',t);
save(filename);


figure(29);
plot(x,y, 'r--')

obstacle = [0,35,8,5];
obstacle2 = [10,10,7,6];
obstacle3 = [10,0,4,8];
obstacle4 = [20,5,7,8];
obstacle5 = [19,14,7,6];
obstacle6 = [8,17,9,7];
obstacle7 = [10,36,7,9];
obstacle8 = [19,30,8,5];
obstacle9 = [20,22,12,4];
obstacle10 = [30,6,5,6];
obstacle11 = [40,6,5,6];
obstacle12 = [38,20,5,6];
obstacle13 = [36,30,5,6];
obstacle14 = [36,30,5,6];
obstacle15 = [36,30,5,6];

initial = [trajectory(end,1) trajectory(end,2)];
target = [trajectory(1,1) trajectory(1,2)];
rectangle('Position',obstacle,'FaceColor',[0 0 0])
rectangle('Position',obstacle2,'FaceColor',[0 0 0])
rectangle('Position',obstacle3,'FaceColor',[0 0 0])
rectangle('Position',obstacle5,'FaceColor',[0 0 0])
rectangle('Position',obstacle6,'FaceColor',[0 0 0])
rectangle('Position',obstacle7,'FaceColor',[0 0 0])
rectangle('Position',obstacle8,'FaceColor',[0 0 0])
rectangle('Position',obstacle9,'FaceColor',[0 0 0])
rectangle('Position',obstacle12,'FaceColor',[0 0 0])
rectangle('Position',obstacle13,'FaceColor',[0 0 0])
rectangle('Position',obstacle14,'FaceColor',[0 0 0])
rectangle('Position',obstacle15,'FaceColor',[0 0 0])

plot(trajectory(end,1),trajectory(end,2),'ko','MarkerFaceColor','g','MarkerSize',8);
plot(trajectory(:,1),trajectory(:,2),'k.');

r = 2; n = 100;
[X,Y,Z] = cylinder(r,n);
h = -0.1;
X = X+trajectory(1,1);
Y = Y+trajectory(1,2);
Z = Z*h;
fill3(X(2,:),Y(2,:),Z(2,: ...
    ),[200/255 1 1],'EdgeColor','k')
