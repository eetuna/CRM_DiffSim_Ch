theta_input = angles(1,3:data_length)';
tau_v1 = tau_v(1,:)';
z = [tau_v1,theta_input ];

z = iddata(tau_v1,theta_input,delta_T,'Name','Dynamics','TimeUnit','s');
sys_1 = arx(z, [na, nb, nk], 'Ts', delta_T );
% compare(z, sys_1);
y1_theta = [sys_1.A, sys_1.B];


phi_input = angles(2,3:data_length)';
tau_v1 = tau_v(1,:)';
z = [tau_v1,phi_input ];

z = iddata(tau_v1,phi_input,delta_T,'Name','Dynamics','TimeUnit','s');
sys_2 = arx(z, [na, nb, nk], 'Ts', delta_T );

y1_phi = [sys_2.A, sys_2.B];


theta_dot_input = anguler_velocties(1,3:data_length-1)';
tau_v1 = tau_v(1,:)';
z = [tau_v1,theta_dot_input ];

z = iddata(tau_v1,theta_dot_input,delta_T,'Name','Dynamics','TimeUnit','s');
sys_3 = arx(z, [na, nb, nk], 'Ts', delta_T );
% compare(z, sys_1);
y1_thetadot = [sys_3.A, sys_3.B];


phi_dot_input = anguler_velocties(2,3:data_length-1)';
tau_v1 = tau_v(1,:)';
z = [tau_v1, phi_dot_input ];

z = iddata(tau_v1,phi_dot_input,delta_T,'Name','Dynamics','TimeUnit','s');
sys_4 = arx(z, [na, nb, nk], 'Ts', delta_T );

y1_phidot = [sys_4.A, sys_4.B];