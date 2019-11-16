function [J_b, Rank] = body_jac(theta)    
%A function to calculate the body jacobian of the intellidex robot.
    %Segment lengths, in inches
    L0 = -13; L1 = 14.7; L2 = 12; L3 = 12; L4 = 9;

    %points along the manipulator
    q_sh = [0 L0 L1].'; %"shoulder"
    q_e  = [L2 L0 L1].'; %"elbow"
    q_w  = [L2+L3 L0 L1].'; %"wrist"
    q_t  = [L2+L3+L4 L0 L1].'; %"tool"

    q1 = q_sh; q2 = q_sh; q3 = q_sh;
    q4 = q_e;
    q5 = q_w;
    q6 = q_t;

    %array of omegas
    w1 = [0 0 1].';
    w2 = [0 -1 0].';
    w3 = [0 0 -1].';
    w4 = [0 0 -1].';
    w5 = [0 0 -1].';
    w6 = [1 0 0].';

    %Calculate the 3-by-3 rotation matrices
    ew1 = Rodriguez(w1, theta(1));
    ew2 = Rodriguez(w2, theta(2));
    ew3 = Rodriguez(w3, theta(3));
    ew4 = Rodriguez(w4, theta(4));
    ew5 = Rodriguez(w5, theta(5));
    ew6 = Rodriguez(w6, theta(6));

    %Calculate the 4-by-4 translation matrices
    ep1 = [ew1 q1; 0 0 0 1];
    ep2 = [ew2 q2; 0 0 0 1];
    ep3 = [ew3 q3; 0 0 0 1];
    ep4 = [ew4 q4; 0 0 0 1];
    ep5 = [ew5 q5; 0 0 0 1];
    ep6 = [ew6 q6; 0 0 0 1];

    %zero configuration transformation, gst(0)
    gst_0 = [eye(3) q_t; 0 0 0 1];

    %Calculating the 6-by-1 twist coordinates
    Psi_1 = [cross(w1,q1); w1];
    Psi_2 = [cross(w2,q2); w2];
    Psi_3 = [cross(w3,q3); w3];
    Psi_4 = [cross(w4,q4); w4];
    Psi_5 = [cross(w5,q5); w5];
    Psi_6 = [cross(w6,q6); w6];

    %calculating the columns of the body jacobian
    Psi_1d = adjoint_inv(ep1*ep2*ep3*ep4*ep5*ep6*gst_0)*Psi_1;
    Psi_2d = adjoint_inv(ep2*ep3*ep4*ep5*ep6*gst_0)*Psi_2;
    Psi_3d = adjoint_inv(ep3*ep4*ep5*ep6*gst_0)*Psi_3;
    Psi_4d = adjoint_inv(ep4*ep5*ep6*gst_0)*Psi_4;
    Psi_5d = adjoint_inv(ep5*ep6*gst_0)*Psi_5;
    Psi_6d = adjoint_inv(ep6*gst_0)*Psi_6;

    J_b = [Psi_1d Psi_2d Psi_3d Psi_4d Psi_5d Psi_6d];

    tol = 0.01; %tolerance when determining rank of J_b
    Rank = rank(J_b, tol);
% 
%     if Rank < 6
%         fprintf("This configuration is singular! Rank = %i\n", Rank);
%     end
end