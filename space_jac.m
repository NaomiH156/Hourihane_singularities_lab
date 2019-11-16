function [J_s, Rank] = space_jac(theta)
%A function to calculate the spatial jacobian of the intellidex robot.
    %Segment lengths, in inches
    L0 = -13; L1 = 14.7; L2 = 12; L3 = 12; L4 = 9;

    %points along the manipulator
    q_sh = [0 L0 L1].'; %"shoulder"
    q_e  = [L2 L0 L1].'; %"elbow"
    q_w  = [L2+L3 L0 L1].'; %"wrist"
    q_t  = [L2+L3+L4 L0 L1].'; %"tool"

    %translating those into an indexed array
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

    %theta = [0 pi/3 0 pi/4 pi/3 pi/12];
    %theta = [pi/3 0 -pi/4 pi/3 0 pi/6];

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


    %calculate the "omega-prime"s
    w2p = ew1*w2;
    w3p = ew1*ew2*w3;
    w4p = ew1*ew2*ew3*w4;
    w5p = ew1*ew2*ew3*ew4*w5;
    w6p = ew1*ew2*ew3*ew4*ew5*w6;

    %calculate the "q-prime"s
    q2p = [q1; 1] + ep1*[q2; 1];
    q3p = q2p + ep1*ep2*[q3; 1];
    q4p = q3p + ep1*ep2*ep3*[q4; 1];
    q5p = q4p + ep1*ep2*ep3*ep4*[q5; 1];
    q6p = q5p + ep1*ep2*ep3*ep4*ep5*[q6; 1];

    %calculate the columns of the body jacobian
    Psi_1 = [cross(w1, q1); w1];
    Psi_2p = [cross(w2p, q2p(1:3)); w2p];
    Psi_3p = [cross(w3p, q3p(1:3)); w3p];
    Psi_4p = [cross(w4p, q4p(1:3)); w4p];
    Psi_5p = [cross(w5p, q5p(1:3)); w5p];
    Psi_6p = [cross(w6p, q6p(1:3)); w6p];

    J_s = [Psi_1 Psi_2p Psi_3p Psi_4p Psi_5p Psi_6p];

    tol = 0.01; %tolerance when determining rank of J_s
    Rank = rank(J_s, tol);
% 
%     if Rank < 6
%         fprintf("This configuration is singular! Rank = %i\n", Rank);
%     end
end
