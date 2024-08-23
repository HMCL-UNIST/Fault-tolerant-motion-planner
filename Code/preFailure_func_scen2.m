%% Pre-failure
function [PFTC_state,PFTC_state_dot,PFTC_thrust,pre_err_norm,pre_err_norm_pos,time_sim,...
          ax25,ax25_enlarge,ax26,ax26_enlarge,ax27,ax27_enlarge] = preFailure_func
    sim('PreFailure');

    PFTC_state = PFTC_data.signals.values;
    PFTC_state_dot = PFTC_data_dot.signals.values;
    PFTC_thrust = PFTC_force.signals.values;

    PFTC_x = PFTC_state(:,1);    PFTC_y = PFTC_state(:,2);      PFTC_z = PFTC_state(:,3);
    PFTC_phi = PFTC_state(:,10); PFTC_theta = PFTC_state(:,11);

    phi_s = deg2rad(0);
    theta_s = 0;

    for i = 1:length(PFTC_x)
        pre_err_norm(i) = 0;
    end
    for i = 1:length(PFTC_x)
        pre_err_norm_pos(i) = norm([PFTC_x(i)-10, PFTC_y(i)-0, PFTC_z(i)-3]);
    end

    figure(25); hold on;
    plot(time_sim,PFTC_thrust(:,1),'k--','LineWidth',1)
    plot(time_sim, PFTC_thrust(:,2), 'k:','LineWidth',1)
    ax25= gca;

    axes('Position',[.6 .32 .2 .25]);
    box on; hold on; grid on;
    plot(time_sim, PFTC_thrust(:,2), 'k:','LineWidth',1)
    ax25_enlarge = gca;

    figure(26); hold on;
    plot(time_sim,PFTC_thrust(:,1),'k--','LineWidth',1)
    plot(time_sim, PFTC_thrust(:,2), 'k:','LineWidth',1)
    ax26= gca;

    axes('Position',[.6 .32 .2 .25]);
    box on; hold on; grid on;
    plot(time_sim, PFTC_thrust(:,2), 'k:','LineWidth',1)
    ax26_enlarge = gca;

    figure(27); hold on;
    plot(time_sim,PFTC_thrust(:,1),'k--','LineWidth',1)
    plot(time_sim, PFTC_thrust(:,2), 'k:','LineWidth',1)
    ax27 = gca;

    axes('Position',[.6 .32 .2 .25]);
    box on; hold on; grid on;
    plot(time_sim, PFTC_thrust(:,2), 'k:','LineWidth',1)
    ax27_enlarge = gca;
    
    figure(28); hold on; grid on
    plot3(0,0,3,'ko','MarkerFaceColor','g','MarkerSize',10);
    plot3(PFTC_x,PFTC_y,PFTC_z,'k-','LineWidth',1)
    PFTC_thrust(:,1);
end