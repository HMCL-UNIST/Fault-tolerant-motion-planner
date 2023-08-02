%% Pre-failure
function [PFTC_state,PFTC_state_dot,PFTC_thrust,pre_err_norm,pre_err_norm_pos,time_sim] = preFailure_func
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
end