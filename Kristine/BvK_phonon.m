%Born - von Karman phonon calculation
% FCC Cu

Data_100_q_t =   [0.15 0.2  0.25 0.275 0.3  0.35 0.4  0.45 0.5  0.55 0.6  0.65 0.7  0.75 0.8  0.9  1.0]; %measured q in rlu
Data_100_nu_t =  [1.17 1.56 1.92 2.12  2.30 2.64 3.01 3.30 3.62 3.88 4.15 4.34 4.54 4.73 4.86 5.02 5.08]; %measured frequency in THz
Data_100_err_t = [0.04 0.04 0.04 0.04  0.04 0.04 0.04 0.05 0.04 0.05 0.05 0.05 0.05 0.06 0.07 0.07 0.08];
Data_100_q_l =   [0.15 0.2  0.25 0.3  0.4  0.5  0.6  0.65 0.7  0.75 0.8  0.85 0.9  1.0]; %measured q in rlu
Data_100_nu_l =  [1.90 2.42 3.02 3.56 4.47 5.32 6.05 6.37 6.60 6.77 6.99 7.14 7.17 7.19]; %measured frequency in THz
Data_100_err_l = [0.09 0.07 0.08 0.06 0.07 0.07 0.08 0.08 0.08 0.12 0.13 0.14 0.12 0.12];

Data_110_q_t2 =   [0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.9 1.0]; %measured q in rlu
Data_110_nu_t2 =  [1.11 1.69 2.27 2.82 3.37 4.30 5.07 5.71 6.04 6.31 6.54 6.80 7.13 7.19]; %measured frequency in THz
Data_110_err_t2 = [0.03 0.04 0.04 0.04 0.04 0.05 0.06 0.06 0.06 0.07 0.10 0.11 0.15 0.12];
Data_110_q_t1 =   [0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9 1.0]; %measured q in rlu
Data_110_nu_t1 =  [1.35 2.03 2.70 3.34 3.89 4.34 4.55 4.75 5.03 5.08]; %measured frequency in THz
Data_110_err_t1 = [0.04 0.04 0.04 0.04 0.05 0.05 0.05 0.07 0.08 0.08];
Data_110_q_l =    [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9 1.0]; %measured q in rlu
Data_110_nu_l =   [2.03 3.70 5.11 5.97 6.36 6.38 5.91 5.73 5.51 5.19 5.08]; %measured frequency in THz
Data_110_err_l =  [0.10 0.08 0.07 0.08 0.10 0.12 0.10 0.08 0.07 0.07 0.08];

Data_111_q_l =   [0.05 0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5]; %measured q in rlu
Data_111_nu_l =  [1.24 1.86 2.46 2.99 3.59 4.54 5.43 6.14 6.67 7.06 7.25 7.40]; %measured frequency in THz
Data_111_err_l = [0.06 0.07 0.07 0.06 0.06 0.06 0.07 0.07 0.08 0.10 0.13 0.13];
Data_111_q_t =   [0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5]; %measured q in rlu
Data_111_nu_t =  [0.79 1.01 1.23 1.47 1.87 2.29 2.66 2.97 3.17 3.34 3.37]; %measured frequency in THz
Data_111_err_t = [0.04 0.05 0.06 0.06 0.06 0.05 0.06 0.06 0.07 0.07 0.07];

Data_1k0_q_PI =   [0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0]; %measured q in rlu
Data_1k0_nu_PI =  [7.19 7.17 7.07 6.80 6.44 6.10 5.77 5.47 5.27 5.13 5.08]; %measured frequency in THz
Data_1k0_err_PI = [0.12 0.15 0.14 0.09 0.09 0.08 0.08 0.08 0.08 0.08 0.08];
Data_1k0_q_L =    [0.0  0.1  0.2  0.3  0.4  0.5]; %measured q in rlu
Data_1k0_nu_L =   [5.08 5.03 4.99 4.97 4.89 4.89]; %measured frequency in THz
Data_1k0_err_L =  [0.08 0.08 0.07 0.08 0.09 0.08];

Da2kg = 1.6605e-27; %Dalton to kg converter
M = 63.546*Da2kg; %atomic mass of Cu
dyn2N = 1e-3; %Convert dyn/cm to N/m
a=3.597;    % Cu lattice constant in AA

K_l1 = 2.802e+4*dyn2N;      % Force constant longitudinal nn1
K_t1 = -0.157e+4*dyn2N;   % Force constant transverse nn1
K_l2 = 0.028e+4*dyn2N;      % Force constant longitudinal nn2
K_t2 = -0.017e+4*dyn2N;   % Force constant transverse nn2
K_l3 = 0.111e+4*dyn2N;      % Force constant longitudinal nn3
K_t3 = 0.005e+4*dyn2N;   % Force constant transverse nn3

K_l4 = 0.053e+4*dyn2N;      % Force constant longitudinal nn4
K_t4 = -0.030e+4*dyn2N;   % Force constant transverse nn4


%%%%%%%%%%%%%%%%%%% 1st neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_j = a/2*[1 1 0; 1 0 1; 0 1 1; 1 -1 0; 1 0 -1; 0 1 -1];

Phi_nn1x=[K_l1 0 0; 0 K_t1 0; 0 0 K_t1];
Phi_nn1y=[K_t1 0 0; 0 K_l1 0; 0 0 K_t1];
si = sin(pi/4);
co = cos(pi/4);
Rot1 = [co -si 0; si co 0; 0 0 1];
Phi_1 = Rot1*Phi_nn1x*Rot1^(-1)
Rot2 = [co 0 -si; 0 1 0; si 0 co];
Phi_2 = Rot2*Phi_nn1x*Rot2^(-1);
Rot3 = [1 0 0; 0 co -si; 0 si co];
Phi_3 = Rot3*Phi_nn1y*Rot3^(-1);
Rot4 = [co si 0; -si co 0; 0 0 1];
Phi_4 = Rot4*Phi_nn1x*Rot4^(-1);
Rot5 = [co 0 si; 0 1 0; -si 0 co];
Phi_5 = Rot5*Phi_nn1x*Rot5^(-1);
Rot6 = [1 0 0; 0 co si; 0 -si co];
Phi_6 = Rot6*Phi_nn1y*Rot6^(-1);

%%%%%%%%%%%%%%%%%%% 2nd neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_j = [r_j; a 0 0; 0 a 0; 0 0 a];

Phi_nn2x=[K_l2 0 0; 0 K_t2 0; 0 0 K_t2];
Phi_nn2y=[K_t2 0 0; 0 K_l2 0; 0 0 K_t2];
Phi_nn2z=[K_t2 0 0; 0 K_t2 0; 0 0 K_l2];

Phi_7 = Phi_nn2x;
Phi_8 = Phi_nn2y;
Phi_9 = Phi_nn2z;

%%%%%%%%%%%%%%%%%%% 3rd neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_add = a/2*[2 1 1; 2 -1 -1; 2 -1 1; 2 1 -1; 1 2 1; 1 2 -1; 1 -2 1; 1 -2 -1; 1 1 2; 1 1 -2; 1 -1 2; 1 -1 -2];
r_j = [r_j; r_add];

Phi_nn3x=[K_l3 0 0; 0 K_t3 0; 0 0 K_t3];
Phi_nn3y=[K_t3 0 0; 0 K_l3 0; 0 0 K_t3];

c2 = 2/sqrt(6);
c1 = 1/sqrt(6);
c0 = 1/sqrt(2);
c3 = 1/sqrt(3);
Rot_10 = [c2 c1 c1; 0 c0 -c0; -c3 c3 c3];
Phi_10 = Rot_10^(-1)*Phi_nn3x*Rot_10;
Rot_11 = [c2 -c1 -c1; 0 c0 -c0; c3 c3 c3];
Phi_11 = Rot_11^(-1)*Phi_nn3x*Rot_11;
Rot_12 = [c2 -c1 c1; 0 c0 c0; -c3 -c3 c3];
Phi_12 = Rot_12^(-1)*Phi_nn3x*Rot_12;
Rot_13 = [c2 c1 -c1; 0 c0 c0; c3 -c3 c3];
Phi_13 = Rot_13^(-1)*Phi_nn3x*Rot_13;
Rot_14 = [c1 c2 c1; -c0 0 c0; c3 -c3 c3];
Phi_14 = Rot_14^(-1)*Phi_nn3x*Rot_14;
Rot_15 = [c1 c2 -c1; -c0 0 -c0; -c3 c3 c3];
Phi_15 = Rot_15^(-1)*Phi_nn3x*Rot_15;
Rot_16 = [c1 -c2 c1; -c0 0 c0; -c3 -c3 -c3];
Phi_16 = Rot_16^(-1)*Phi_nn3x*Rot_16;
Rot_17 = [c1 -c2 -c1; -c0 0 -c0; c3 c3 -c3];
Phi_17 = Rot_17^(-1)*Phi_nn3x*Rot_17;
Rot_18 = [c1 c1 c2; c0 -c0 0; c3 c3 -c3];
Phi_18 = Rot_18^(-1)*Phi_nn3x*Rot_18;
Rot_19 = [c1 c1 -c2; -c0 c0 0; c3 c3 c3];
Phi_19 = Rot_19^(-1)*Phi_nn3x*Rot_19;
Rot_20 = [c1 -c1 c2; c0 c0 0; -c3 c3 c3];
Phi_20 = Rot_20^(-1)*Phi_nn3x*Rot_20;
Rot_21 = [c1 -c1 -c2; -c0 -c0 0; -c3 c3 -c3];
Phi_21 = Rot_21^(-1)*Phi_nn3x*Rot_21;

%Phi_nn3_211 = [0.0758 0.0353 0.0353; 0.0353 0.0229 0.0176; 0.0353 0.0176 0.0229]; 

%%%%%%%%%%%%%%%%%%%% 4th neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_add = a/2*[2 2 0; 2 0 2; 0 2 2; 2 -2 0; 2 0 -2; 0 2 -2];
r_j = [r_j; r_add];

Phi_nn4x=[K_l4 0 0; 0 K_t4 0; 0 0 K_t4];
Phi_nn4y=[K_t4 0 0; 0 K_l4 0; 0 0 K_t4];

Phi_22 = Rot1*Phi_nn2x*Rot1^(-1);
Phi_23 = Rot2*Phi_nn2x*Rot2^(-1);
Phi_24 = Rot3*Phi_nn2y*Rot3^(-1);
Phi_25 = Rot4*Phi_nn2x*Rot4^(-1);
Phi_26 = Rot5*Phi_nn2x*Rot5^(-1);
Phi_27 = Rot6*Phi_nn2y*Rot6^(-1);

%%%%%%%%%%%% Calculate dispersion 100 %%%%%%%%%%%%%%%%%%%%%%%%

qstep=0.02;
eq = [1 0 0];   
leq = sqrt(eq*eq');
eq = eq/leq; %unit vertor in q direction
omega_list = [];
q_list = [];
cos_list = [];


for q=-1:qstep:1
    qvec = leq*2*pi/a*q*eq;

    phase = 1-cos(r_j*qvec');
    
    % The next should be a loop over neighbours if I only could ...

    D = 0;
    D = D + Phi_1*phase(1) + Phi_2*phase(2) + Phi_3*phase(3);
    D = D + Phi_4*phase(4) + Phi_5*phase(5) + Phi_6*phase(6);
    D = D + Phi_7*phase(7) + Phi_8*phase(8) + Phi_9*phase(9);
    D = D + Phi_10*phase(10) + Phi_11*phase(11) + Phi_12*phase(12);
    D = D + Phi_13*phase(13) + Phi_14*phase(14) + Phi_15*phase(15);
    D = D + Phi_16*phase(16) + Phi_17*phase(17) + Phi_18*phase(18);
    D = D + Phi_19*phase(19) + Phi_20*phase(20) + Phi_21*phase(21);
    D = D + Phi_22*phase(22) + Phi_23*phase(23) + Phi_24*phase(24);
    D = D + Phi_25*phase(25) + Phi_26*phase(26) + Phi_27*phase(27);

    D = 2*D; %Correct for that we took only half of the neighbours

   [V,DIAG] = eig(D);
   E = diag(DIAG);
  
   omega = sqrt(E/M)/(2*pi*1e12);
   cosines = V'*eq';
   
   om_cos = [omega cosines];
   om_cos_sort = sortrows(om_cos,1,'ascend');
   
   q_list = [q_list q];
   omega_list = [omega_list om_cos_sort(:,1)];
   cos_list = [cos_list abs(om_cos_sort(:,2))];
end

figure(1)
plot(q_list,omega_list(1,:),'r-')
xlabel('q = (h h h) [rlu]')
ylabel('\nu [THz] / direction cosines')
%title('Phonon dispersion in Cu. Data + model: [Svensson67]')
hold on
plot(q_list,omega_list(2,:),'r-')
plot(q_list,omega_list(3,:),'b-')
plot(q_list,cos_list(1,:),'r--')
plot(q_list,cos_list(2,:),'r--')
plot(q_list,cos_list(3,:),'b--')

errorbar(Data_100_q_t, Data_100_nu_t, Data_100_err_t,'rx')
errorbar(Data_100_q_l, Data_100_nu_l, Data_100_err_l,'bo')
hold off
print('-depsc','Cu_phonon_100.eps')

%%%%%%%%%%%% Calculate dispersion 110 %%%%%%%%%%%%%%%%%%%%%%%%

qstep=0.02;
eq = [1 0 1];   
leq = sqrt(eq*eq');
eq = eq/leq; %unit vertor in q direction
omega_list = [];
q_list = [];
cos_list = [];


for q=-1:qstep:1
    qvec = leq*2*pi/a*q*eq;

    phase = 1-cos(r_j*qvec');
    
    % The next should be a loop over neighbours if I only could ...

    D = 0;
    D = D + Phi_1*phase(1) + Phi_2*phase(2) + Phi_3*phase(3);
    D = D + Phi_4*phase(4) + Phi_5*phase(5) + Phi_6*phase(6);
    D = D + Phi_7*phase(7) + Phi_8*phase(8) + Phi_9*phase(9);
    D = D + Phi_10*phase(10) + Phi_11*phase(11) + Phi_12*phase(12);
    D = D + Phi_13*phase(13) + Phi_14*phase(14) + Phi_15*phase(15);
    D = D + Phi_16*phase(16) + Phi_17*phase(17) + Phi_18*phase(18);
    D = D + Phi_19*phase(19) + Phi_20*phase(20) + Phi_21*phase(21);
    D = D + Phi_22*phase(22) + Phi_23*phase(23) + Phi_24*phase(24);
    D = D + Phi_25*phase(25) + Phi_26*phase(26) + Phi_27*phase(27);

    D = 2*D; %Correct for that we took only half of the neighbours

   [V,DIAG] = eig(D);
   E = diag(DIAG);
  
   omega = sqrt(E/M)/(2*pi*1e12);
   cosines = V'*eq';
   
   om_cos = [omega cosines];
   om_cos_sort = sortrows(om_cos,1,'ascend');
   
   q_list = [q_list q];
   omega_list = [omega_list om_cos_sort(:,1)];
   cos_list = [cos_list abs(om_cos_sort(:,2))];
end

figure(1)
plot(q_list,omega_list(1,:),'r-')
xlabel('q = (h h 0) [rlu]')
ylabel('\nu [THz] / direction cosines')
%title('Phonon dispersion in Cu. Data + model: [Svensson67]')
hold on
plot(q_list,omega_list(2,:),'g-')
plot(q_list,omega_list(3,:),'b-')
plot(q_list,cos_list(1,:),'r--')
plot(q_list,cos_list(2,:),'g--')
plot(q_list,cos_list(3,:),'b--')

errorbar(Data_110_q_t2, Data_110_nu_t2, Data_110_err_t2,'gx')
errorbar(Data_110_q_t1, Data_110_nu_t1, Data_110_err_t1,'rx')
errorbar(Data_110_q_l, Data_110_nu_l, Data_110_err_l,'bo')
hold off
print('-depsc','Cu_phonon_110.eps')

%%%%%%%%%%%% Calculate dispersion 111 %%%%%%%%%%%%%%%%%%%%%%%%

qstep=0.02;
eq = [1 1 1];   
leq = sqrt(eq*eq');
eq = eq/leq; %unit vertor in q direction
omega_list = [];
q_list = [];
cos_list = [];


for q=-1:qstep:1
    qvec = leq*2*pi/a*q*eq;

    phase = 1-cos(r_j*qvec');
    
    % The next should be a loop over neighbours if I only could ...

    D = 0;
    D = D + Phi_1*phase(1) + Phi_2*phase(2) + Phi_3*phase(3);
    D = D + Phi_4*phase(4) + Phi_5*phase(5) + Phi_6*phase(6);
    D = D + Phi_7*phase(7) + Phi_8*phase(8) + Phi_9*phase(9);
    D = D + Phi_10*phase(10) + Phi_11*phase(11) + Phi_12*phase(12);
    D = D + Phi_13*phase(13) + Phi_14*phase(14) + Phi_15*phase(15);
    D = D + Phi_16*phase(16) + Phi_17*phase(17) + Phi_18*phase(18);
    D = D + Phi_19*phase(19) + Phi_20*phase(20) + Phi_21*phase(21);
    D = D + Phi_22*phase(22) + Phi_23*phase(23) + Phi_24*phase(24);
    D = D + Phi_25*phase(25) + Phi_26*phase(26) + Phi_27*phase(27);

    D = 2*D; %Correct for that we took only half of the neighbours

   [V,DIAG] = eig(D);
   E = diag(DIAG);
  
   omega = sqrt(E/M)/(2*pi*1e12);
   cosines = V'*eq';
   
   om_cos = [omega cosines];
   om_cos_sort = sortrows(om_cos,1,'ascend');
   
   q_list = [q_list q];
   omega_list = [omega_list om_cos_sort(:,1)];
   cos_list = [cos_list abs(om_cos_sort(:,2))];
end

figure(1)
plot(q_list,omega_list(1,:),'r-')
xlabel('q = (h h h) [rlu]')
ylabel('\nu [THz] / direction cosines')
%title('Phonon dispersion in Cu. Data + model: [Svensson67]')
hold on
plot(q_list,omega_list(2,:),'r-')
plot(q_list,omega_list(3,:),'b-')
plot(q_list,cos_list(1,:),'r--')
plot(q_list,cos_list(2,:),'r--')
plot(q_list,cos_list(3,:),'b--')

errorbar(Data_111_q_t, Data_111_nu_t, Data_111_err_t,'rx')
errorbar(Data_111_q_l, Data_111_nu_l, Data_111_err_l,'bo')
hold off
print('-depsc','Cu_phonon_111.eps')

%%%%%%%%%%%% Calculate dispersion 1k0 %%%%%%%%%%%%%%%%%%%%%%%%

qstep=0.02;

omega_list = [];
q_list = [];
cos_list = [];


for q=-1:qstep:1
    qdir = [1 q 0];
    eq = qdir/sqrt(qdir*qdir'); %unit vector in q direction
    qvec = 2*pi/a*qdir;

    phase = 1-cos(r_j*qvec');
    
    % The next should be a loop over neighbours if I only could ...

    D = 0;
    D = D + Phi_1*phase(1) + Phi_2*phase(2) + Phi_3*phase(3);
    D = D + Phi_4*phase(4) + Phi_5*phase(5) + Phi_6*phase(6);
    D = D + Phi_7*phase(7) + Phi_8*phase(8) + Phi_9*phase(9);
    D = D + Phi_10*phase(10) + Phi_11*phase(11) + Phi_12*phase(12);
    D = D + Phi_13*phase(13) + Phi_14*phase(14) + Phi_15*phase(15);
    D = D + Phi_16*phase(16) + Phi_17*phase(17) + Phi_18*phase(18);
    D = D + Phi_19*phase(19) + Phi_20*phase(20) + Phi_21*phase(21);
    D = D + Phi_22*phase(22) + Phi_23*phase(23) + Phi_24*phase(24);
    D = D + Phi_25*phase(25) + Phi_26*phase(26) + Phi_27*phase(27);

    D = 2*D; %Correct for that we took only half of the neighbours

   [V,DIAG] = eig(D);
   E = diag(DIAG);
  
   omega = sqrt(E/M)/(2*pi*1e12);
   cosines = V'*eq';
   
   om_cos = [omega cosines];
   om_cos_sort = sortrows(om_cos,1,'ascend');
   
   q_list = [q_list q];
   omega_list = [omega_list om_cos_sort(:,1)];
   cos_list = [cos_list abs(om_cos_sort(:,2))];
end

figure(1)
plot(q_list,omega_list(1,:),'r-')
xlabel('q = (1 k 0) [rlu]')
ylabel('\nu [THz] / direction cosines')
%title('Phonon dispersion in Cu. Data + model: [Svensson67]')
hold on
plot(q_list,omega_list(2,:),'r-')
plot(q_list,omega_list(3,:),'b-')
plot(q_list,cos_list(1,:),'r--')
plot(q_list,cos_list(2,:),'r--')
plot(q_list,cos_list(3,:),'b--')

errorbar(Data_1k0_q_PI, Data_1k0_nu_PI, Data_1k0_err_PI,'rx')
errorbar(Data_1k0_q_L, Data_1k0_nu_L, Data_1k0_err_L,'bo')

print('-depsc','Cu_phonon_1k0.eps')
hold off