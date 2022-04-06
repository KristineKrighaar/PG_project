%Born - von Karman phonon calculation
% Hexagonal Pyrolytic Graphite

q2q = 0.1/114;  % Conversion from mm on zoomed picture to rlu
E2E = 2/91;     % Conversion from mm on zoomed picture to THz

Data_100_q_ta =       [0.02 0.04 0.06 0.08 0.10 0.12 0.145 0.20 0.25 0.295]; %measured q in rlu
Data_100_nu_ta =  E2E*[6    12   20   28   42   55   80    136  209  294]; %measured frequency in mm->THz
Data_100_err_ta = E2E*[3    3    3    3    3    3    3     3    3    3];
Data_100_q_to =       [0.00 0.05 0.10 0.145 0.20 0.245]; %measured q in rlu
Data_100_nu_to =  E2E*[171  173  177  190   218  264]; %measured frequency in mm->THz
Data_100_err_to = E2E*[3    3    3    3     3    3];
Data_100_q_la =   q2q*[14 13 17 16  28  34  51  57  96]; %measured q in mm->rlu
Data_100_nu_la =  E2E*[66 81 88 135 135 142 270 270 447]; %measured E in mm->THz
Data_100_err_la = E2E*[3  3  3  3   3   3   3   3   3];
Data_100_q_lo =   q2q*[0  10 19  22  48  46  66  86]; % measured q in mm->rlu
Data_100_nu_lo =  E2E*[68 91 112 134 225 270 347 447]; %measured E in mm->THz
Data_100_err_lo = E2E*[3  3  3   3   3   3   3   3];


Data_001_q_la =   [0.1 0.2 0.3 0.4 0.5]; %measured q in rlu
Data_001_nu_la =  E2E*[30 50 75 100 120]; %measured frequency in mm->THz
Data_001_err_la = E2E*[3 3 3 3 3];
Data_001_q_ta =   [0.3 0.4 0.5]; %measured q in rlu
Data_001_nu_ta =  E2E*[29 40 48]; %measured frequency in mm->THz
Data_001_err_ta = E2E*[3 3 3];
Data_001_q_lo =   [0.0 0.1 0.2 0.3 0.4 0.5]; %measured q in rlu
Data_001_nu_lo =  E2E*[171 172 166 152 138 120]; %measured frequency in mm->THz
Data_001_err_lo = E2E*[3 3 3 3 3 3];
Data_001_q_to =   [0.0 0.1 0.2 0.3 0.4 0.5]; %measured q in rlu
Data_001_nu_to =  E2E*[68 68 65 62 54 48]; %measured frequency in mm->THz
Data_001_err_to = E2E*[3 3 3 3 3 3];

Da2kg = 1.6605e-27; %Dalton to kg converter
M = 12.011*Da2kg; %atomic mass of C
dyn2N = 1e-3; %Convert dyn/cm to N/m
a=2.461;    % PG lattice constant in AA  paper gives 2.45
c=6.708;    % PG lattice constant in AA  paper gives 6.70

K_l1 =  3.62e+5*dyn2N;     % Force constant longitudinal nn1 - in (a,b) plane
K_t1 =  1.99e+5*dyn2N;     % Force constant transverse nn1 - in (a,b) plane
K_l2 =  1.33e+5*dyn2N;     % Force constant longitudinal nn2 - in (a,b) plane
K_t2 = -0.520e+5*dyn2N;    % Force constant transverse nn2 - in (a,b) plane
K_l3 = -0.037e+5*dyn2N;    % Force constant longitudinal nn3 - in (a,b) plane
K_t3 =  0.288e+5*dyn2N;    % Force constant transverse nn3 - in (a,b) plane

K_l4 =  0.058e+5*dyn2N;    % Force constant longitudinal nn4 - along c
K_t4 =  0.0077e+5*dyn2N;   % Force constant transverse nn4 - along c

%%%%%%%%%%%%%%%%%%%%%%%%% Lattice basis %%%%%%%%%%%%%%%%%%%%%%%%%%%

% PG is hexagonal close packed with 4 atoms per cell
% Coordinates from the paper Fig. 7: Atoms A, B, C, D
s3 = sqrt(3);
Delta = [0 0 0; a/(2*s3) a/2 0; 0 0 c/2; -a/(2*s3) a/2 c/2];

% Lattice vectors from paper fig. 7
avec = a*[s3/2 0.5 0];
bvec = a*[-s3/2 0.5 0];
cvec = c*[0 0 1];

% Rotation matrices
Rot120 = [-0.5 s3/2 0; -s3/2 -0.5 0; 0 0 1];
Rot60  = [0.5 s3/2 0; -s3/2 0.5 0; 0 0 1];
%%%%%%%%%%%%%%%%%%% 1st neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


r_j1 = a*[-1/s3 0 0; 0.5/s3 -0.5 0; 0.5/s3 0.5 0]; % This holds from atoms A and D
r_j2 = -r_j1; % This holds from atoms B and C

Phi_nn1=[K_l1 0 0; 0 K_t1 0; 0 0 K_t1];
Phi1 = Phi_nn1;

Rot2=Rot120;
Phi2 = Rot2*Phi_nn1*Rot2^(-1);
Rot3 = Rot120*Rot120;
Phi3 = Rot3*Phi_nn1*Rot3^(-1);

astar = 4*pi/(sqrt(3)*a);  % length of reciprocal lattice vector a*
cstar = 2*pi/c;  % length of reciprocal lattice vector c*



%%%%%%%%%%%%%%%%%%% 2nd neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_add = [0 a 0; -a*s3/2 a/2 0; -a*s3/2 -a/2 0; 0 -a 0; a*s3/2 -a/2 0; a*s3/2 a/2 0];
r_j1 = [r_j1; r_add];
r_j2 = [r_j2; r_add];
% 
Phi_nn2=[K_t2 0 0; 0 K_l2 0; 0 0 K_t2];
 
Phi4 = Phi_nn2;
Rot5 = Rot60;
Phi5 = Rot5*Phi_nn2*Rot5^(-1);
Rot6 = Rot120;
Phi6 = Rot6*Phi_nn2*Rot6^(-1);
Rot7 = Rot60*Rot120;
Phi7 = Rot7*Phi_nn2*Rot7^(-1);
Rot8 = Rot120*Rot120;
Phi8 = Rot8*Phi_nn2*Rot8^(-1);
Rot9 = Rot60^(-1);
Phi9 = Rot9*Phi_nn2*Rot9^(-1);

%%%%%%%%%%%%%%%%%%% 3rd neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_add = [2*a/s3 0 0; -a/s3 a 0; -a/s3 -a 0];
r_j1 = [r_j1; r_add];
r_j2 = [r_j2; -r_add];

 
Phi_nn3=[K_l3 0 0; 0 K_t3 0; 0 0 K_t3];

% c2 = 2/sqrt(6);
% c1 = 1/sqrt(6);
% c0 = 1/sqrt(2);
% c3 = 1/sqrt(3);

Phi10 = Phi_nn3;
Rot11 = Rot120;
Phi11 = Rot11*Phi_nn3*Rot11^(-1);
Rot12 = Rot120*Rot120;
Phi12 = Rot12*Phi_nn3*Rot12^(-1);


% %%%%%%%%%%%%%%%%%%%% 4th neighbour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_add = c/2*[0 0 1; 0 0 -1];
r_j1 = [r_j1; r_add]; % Sublattice B and D do not have this coupling

Phi_nn4=[K_t4 0 0; 0 K_t4 0; 0 0 K_l4];
 
Phi13 = Phi_nn4;
Phi14 = Phi_nn4;


%%%%%%%%%%%% Calculate dispersion 100 %%%%%%%%%%%%%%%%%%%%%%%%

qstep=0.005;
eq = [1 0 0];   
leq = sqrt(eq*eq');
eq = eq/leq; %unit vector in q direction
omega_list = [];
q_list = [];
cos_list = [];

Zero_3Dim = zeros(3,3);

for q=0:qstep:0.5

    qvec = leq*astar*q*eq;

    
% Make the matrix in the 6D space [r_A r_B]

Phi_diag  = Phi1 +Phi2 +Phi3 +Phi4*(1-exp(i*qvec*r_j1(4,:)'))+Phi5*(1-exp(i*qvec*r_j1(5,:)')) ...
                             +Phi6*(1-exp(i*qvec*r_j1(6,:)'))+Phi7*(1-exp(i*qvec*r_j1(7,:)')) ...
                             +Phi8*(1-exp(i*qvec*r_j1(8,:)'))+Phi9*(1-exp(i*qvec*r_j1(9,:)')) ...
           +Phi10+Phi11+Phi12;
Phi_diag1 = Phi13+Phi14;
Phi_offdiag1 = -Phi1*exp(i*qvec*r_j1(1,:)')  -Phi2*exp(i*qvec*r_j1(2,:)')  -Phi3*exp(i*qvec*r_j1(3,:)') ...
               -Phi10*exp(i*qvec*r_j1(10,:)')-Phi11*exp(i*qvec*r_j1(11,:)')-Phi12*exp(i*qvec*r_j1(12,:)');
Phi_offdiag2 = -Phi1*exp(i*qvec*r_j2(1,:)')  -Phi2*exp(i*qvec*r_j2(2,:)')  -Phi3*exp(i*qvec*r_j2(3,:)') ...
               -Phi10*exp(i*qvec*r_j2(10,:)')-Phi11*exp(i*qvec*r_j2(11,:)')-Phi12*exp(i*qvec*r_j2(12,:)'); %identical the to complex conjugate of Phi_offdiag1

Phi_6Dim = [Phi_diag+Phi_diag1 Phi_offdiag1; Phi_offdiag2 Phi_diag];
    
Phi_AC = -Phi13*exp(i*qvec*r_j1(13,:)') -Phi14*exp(i*qvec*r_j1(14,:)');

Phi_6Dim_offdiag = [Phi_AC Zero_3Dim; Zero_3Dim Zero_3Dim];

% Make the matrix in the 12D space [r_A r_B r_C r_D]

Phi_12Dim = [Phi_6Dim Phi_6Dim_offdiag; conj(Phi_6Dim_offdiag) Phi_6Dim];

   [V,DIAG] = eig(Phi_12Dim);
   E = diag(DIAG);
  
  omega = sqrt(E/M)/(2*pi*1e12);
  cosines = V'*[eq' ; eq'; eq'; eq'];

   om_cos = [omega cosines];
   om_cos_sort = sortrows(om_cos,1,'ascend');
   
   q_list = [q_list q];
   omega_list = [omega_list om_cos_sort(:,1)];
   cos_list = [cos_list abs(om_cos_sort(:,2))];
end

%omega_minus = sqrt(K_t1/M*(3-sqrt(5+4*cos(2*pi*q_list))))/(2*pi*1e12); % z-polarization modes calculated by hand for nn only
%omega_plus = sqrt(K_t1/M*(3+sqrt(5+4*cos(2*pi*q_list))))/(2*pi*1e12);

figure(1)
plot(q_list,omega_list(1,:),'r-')
axis([0 0.5 0 12])
xlabel('q = (h 0 0) [rlu]')
ylabel('\nu [THz] / direction cosines')
title('Phonon dispersion in PG. Data + model: [Nicklow72]')
hold on
plot(q_list,omega_list(2,:),'k-')
plot(q_list,omega_list(3,:),'r-')
plot(q_list,omega_list(4,:),'b-')
plot(q_list,omega_list(5,:),'g-')
plot(q_list,omega_list(6,:),'k-')
plot(q_list,omega_list(7,:),'k-')
plot(q_list,omega_list(8,:),'r-')
plot(q_list,omega_list(9,:),'b-')
plot(q_list,omega_list(10,:),'g-')
plot(q_list,omega_list(11,:),'k-')
plot(q_list,omega_list(12,:),'g-')


%plot(q_list,omega_minus,'k--')
%plot(q_list,omega_plus,'g--')

% plot(q_list,cos_list(1,:),'r--')
% plot(q_list,cos_list(2,:),'r--')
% plot(q_list,cos_list(3,:),'b--')
% 
errorbar(Data_100_q_ta, Data_100_nu_ta, Data_100_err_ta,'rx')
errorbar(Data_100_q_to, Data_100_nu_to, Data_100_err_to,'bx')
errorbar(Data_100_q_la, Data_100_nu_la, Data_100_err_la,'ro')
errorbar(Data_100_q_lo, Data_100_nu_lo, Data_100_err_lo,'bo')
hold off

print('-depsc','PG_phonon_100.eps')


%%%%%%%%%%%% Calculate dispersion 001 %%%%%%%%%%%%%%%%%%%%%%%%

qstep=0.005;
eq = [0 0 1];   
leq = sqrt(eq*eq');
eq = eq/leq; %unit vector in q direction
omega_list = [];
q_list = [];
cos_list = [];

Zero_3Dim = zeros(3,3);

for q=0:qstep:0.5

    qvec = leq*cstar*q*eq;

    
% Make the matrix in the 6D space [r_A r_B]

Phi_diag  = Phi1 +Phi2 +Phi3 +Phi4*(1-exp(i*qvec*r_j1(4,:)'))+Phi5*(1-exp(i*qvec*r_j1(5,:)')) ...
                             +Phi6*(1-exp(i*qvec*r_j1(6,:)'))+Phi7*(1-exp(i*qvec*r_j1(7,:)')) ...
                             +Phi8*(1-exp(i*qvec*r_j1(8,:)'))+Phi9*(1-exp(i*qvec*r_j1(9,:)')) ...
           +Phi10+Phi11+Phi12;
Phi_diag1 = Phi13+Phi14;
Phi_offdiag1 = -Phi1*exp(i*qvec*r_j1(1,:)')  -Phi2*exp(i*qvec*r_j1(2,:)')  -Phi3*exp(i*qvec*r_j1(3,:)') ...
               -Phi10*exp(i*qvec*r_j1(10,:)')-Phi11*exp(i*qvec*r_j1(11,:)')-Phi12*exp(i*qvec*r_j1(12,:)');
Phi_offdiag2 = -Phi1*exp(i*qvec*r_j2(1,:)')  -Phi2*exp(i*qvec*r_j2(2,:)')  -Phi3*exp(i*qvec*r_j2(3,:)') ...
               -Phi10*exp(i*qvec*r_j2(10,:)')-Phi11*exp(i*qvec*r_j2(11,:)')-Phi12*exp(i*qvec*r_j2(12,:)'); %identical the to complex conjugate of Phi_offdiag1

Phi_6Dim = [Phi_diag+Phi_diag1 Phi_offdiag1; Phi_offdiag2 Phi_diag];
    
Phi_AC = -Phi13*exp(i*qvec*r_j1(13,:)') -Phi14*exp(i*qvec*r_j1(14,:)');

Phi_6Dim_offdiag = [Phi_AC Zero_3Dim; Zero_3Dim Zero_3Dim];

% Make the matrix in the 12D space [r_A r_B r_C r_D]

Phi_12Dim = [Phi_6Dim Phi_6Dim_offdiag; conj(Phi_6Dim_offdiag) Phi_6Dim];

   [V,DIAG] = eig(Phi_12Dim);
   E = diag(DIAG);
  
  omega = sqrt(E/M)/(2*pi*1e12);
  cosines = V'*[eq' ; eq'; eq'; eq'];

   om_cos = [omega cosines];
   om_cos_sort = sortrows(om_cos,1,'ascend');
   
   q_list = [q_list q];
   omega_list = [omega_list om_cos_sort(:,1)];
   cos_list = [cos_list abs(om_cos_sort(:,2))];
end

%omega_minus = sqrt(K_t1/M*(3-sqrt(5+4*cos(2*pi*q_list))))/(2*pi*1e12); % z-polarization modes calculated by hand for nn only
%omega_plus = sqrt(K_t1/M*(3+sqrt(5+4*cos(2*pi*q_list))))/(2*pi*1e12);

q_list = -q_list; % Flip x-axis to compare with article

figure(2)
plot(q_list,omega_list(1,:),'r-')
axis([-0.5 0 0 12])
xlabel('q = (0 0 l) [rlu]')
ylabel('\nu [THz] / direction cosines')
title('Phonon dispersion in PG. Data + model: [Nicklow72]')
hold on
plot(q_list,omega_list(2,:),'k-')
plot(q_list,omega_list(3,:),'r-')
plot(q_list,omega_list(4,:),'b-')
plot(q_list,omega_list(5,:),'g-')
plot(q_list,omega_list(6,:),'k-')
plot(q_list,omega_list(7,:),'k-')
plot(q_list,omega_list(8,:),'r-')
plot(q_list,omega_list(9,:),'b-')
plot(q_list,omega_list(10,:),'g-')
plot(q_list,omega_list(11,:),'k-')
plot(q_list,omega_list(12,:),'g-')


% plot(q_list,cos_list(1,:),'r--')
% plot(q_list,cos_list(2,:),'r--')
% plot(q_list,cos_list(3,:),'b--')
% 
errorbar(-Data_001_q_ta, Data_001_nu_ta, Data_001_err_ta,'rx')
errorbar(-Data_001_q_to, Data_001_nu_to, Data_001_err_to,'bx')
errorbar(-Data_001_q_la, Data_001_nu_la, Data_001_err_la,'ro')
errorbar(-Data_001_q_lo, Data_001_nu_lo, Data_001_err_lo,'bo')
hold off

print('-depsc','PG_phonon_001.eps')
