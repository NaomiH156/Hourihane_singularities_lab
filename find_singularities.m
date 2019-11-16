%This script calls functions that calculate the spatial and body jacobians.
clear all

theta1 = [0 pi/3 0 pi/4 pi/3 pi/12];
theta2 = [pi/3 0 -pi/4 pi/3 0 pi/6];

[J_s_1, Rank_space_1] = space_jac(theta1)
[J_b_1, Rank_body_1] = body_jac(theta1)

[J_s_2, Rank_space_2] = space_jac(theta2)
[J_b_2, Rank_body_2] = body_jac(theta2)
%theta 2 is singular because its body and spatial jacobians have a rank of
    %5.

%Determine the combination of joint movements that produce no end-effector
    %motion, and set of end effector velocities that cannot be achieved.
    %The set of joint movements that produce no end-effector velocity is 
        %simply the null space of the Jacobian.
null_body_2 = null(J_b_2)
null_space_2 = null(J_s_2)

%these represent the set of joint velocities that produce little or no end
    %effector motion at the singular configuration. Show that they produce
    %no velocity:
V_space = J_s_2 * null_space_2
V_body = J_b_2 * null_body_2

%These velocities are on the order of 10^-13 inches per second. Thus, at
%this configuration, ANY velocity cannot be achieved. The robot is "stuck."

%Next, pick a singular configuration. Determine the internal motion space
    %and set of end effector velocities that cannot be achieved.
theta3 = [0 pi/2 0 0 0 0];
[J_s_3, Rank_space_3] = space_jac(theta3)
[J_b_3, Rank_body_3] = body_jac(theta3)

%This is a singular configuration because the body jacobian has a rank of 5
null_body_3 = null(J_b_3)
    %These are the joint velocities that produce no end effector motion
V_body_3 = J_b_3*null_body_3
