% =========================================================================
% Script di Inizializzazione e Linearizzazione: Elicottero 2DoF
% Progetto di Controllo Robusto (Sintesi Robusta su 4 stati)
% =========================================================================
clear all; close all; clc;
disp('Configurazione parametri Elicottero 2DoF (Controllo Robusto)...');
which uss -all

%% 4. PUNTO DI EQUILIBRIO (Trim Point)

alpha_0 = 0;
beta_0  = 0;

% Equilibrio orizzontale senza coppie applicate
F1_0 = 0;
F2_0 = 0;

%% 1. PARAMETRI NOMINALI (Perfettamente noti)
p.J_y = 0.00023; 
p.J_z = 0.00364; 
p.I_b = 0.00023; 
p.m = 0.2; 
p.l_nom = 0.2; 
p.c_alpha = 0.01; 
p.c_beta = 0.01;
p.g = 9.81;

%% PARAMETRI NOMINALI PER LQG

J_alpha_nom = 0.012;
l_nom       = 0.2;
eps_p_nom   = 0.1;
eps_y_nom   = 0.1;
tau_d_nom   = 0.02;
Jbeta_nom = p.J_y*sin(alpha_0)^2 + ...
            (p.J_z+p.m*l_nom^2)*cos(alpha_0)^2 + ...
            p.I_b;


%% PARAMETRI INCERTI PER ROBUSTEZZA

J_alpha_u = ureal('J_alpha',0.012,'Percentage',10);
l_u       = ureal('l',0.2,'Percentage',5);
eps_p_u   = ureal('eps_p',0.1,'Percentage',20);
eps_y_u   = ureal('eps_y',0.1,'Percentage',20);
tau_d_u   = ureal('tau_d',0.02,'Percentage',20);
Jy_u = ureal('Jy',0.00023,'Percentage',5);
Jz_u = ureal('Jz',0.00364,'Percentage',5);
m_u  = ureal('m',0.2,'Percentage',5);

%% 3. DINAMICA DEGLI ATTUATORI
omega_n = 30;   
zeta    = 0.7;     

% Vettori per il secondo ordine
num_2nd = omega_n^2;
den_2nd = [1, 2*zeta*omega_n, omega_n^2];


num_pade_nom = [-tau_d_nom/2 1];
den_pade_nom = [tau_d_nom/2 1];

G_delay_pade_nom = tf(num_pade_nom,den_pade_nom);


num_pade_unc = [-tau_d_u/2 1];
den_pade_unc = [tau_d_u/2 1];

G_delay_pade_unc = tf(num_pade_unc,den_pade_unc);

G_act_2nd = tf(num_2nd,den_2nd);
G_act_nom = G_act_2nd*G_delay_pade_nom;
G_actuator_unc = G_act_2nd*G_delay_pade_unc;


%% Stato e ingressi simbolici

syms alpha alphad beta betad F1 F2 real
syms J_alpha_s l_s eps_p_s eps_y_s real
syms Jy_s Jz_s m_s real


x = [alpha;
     alphad;
     beta;
     betad];

u = [F1;
     F2];


%% Inerzia yaw variabile

% valori nominali per linearizzazione

p.J_alpha = 0.012;


Jbeta = Jy_s*sin(alpha)^2 + ...
        (Jz_s+m_s*l_s^2)*cos(alpha)^2 + ...
        p.I_b;

alphadd = (-p.c_alpha*alphad ...
          -p.m*p.g*l_s*sin(alpha) ...
          +l_s*F1*cos(beta) ...
          +eps_p_s*l_s*F2*sin(beta)) / J_alpha_s;


betadd = (-p.c_beta*betad ...
          +l_s*F2*cos(alpha) ...
          +eps_y_s*l_s*F1*sin(alpha)) / Jbeta;


f = [
    alphad;
    alphadd;
    betad;
    betadd
];


%% Linearizzazione

A_sym = jacobian(f,x);
B_sym = jacobian(f,u);




%% Punto di equilibrio

x0 = [
alpha_0;
0;
beta_0;
0
];


u0 = [
F1_0;
F2_0
];

%% LINEARIZZAZIONE NOMINALE

f_nom = subs(f,...
[J_alpha_s,l_s,eps_p_s,eps_y_s,Jy_s,Jz_s,m_s],...
[J_alpha_nom,l_nom,eps_p_nom,eps_y_nom,p.J_y,p.J_z,p.m]);

A_sym_nom = jacobian(f_nom,x);
B_sym_nom = jacobian(f_nom,u);


A_nom = double(subs(A_sym_nom,...
[x;u],...
[x0;u0]));


B_nom = double(subs(B_sym_nom,...
[x;u],...
[x0;u0]));



%% SENSORI IMU - Accelerometro + Magnetometro

h = [
    -p.g*sin(alpha);
    cos(alpha)*cos(beta);
    -sin(beta)
];

C_sym = jacobian(h,x);

Cy = double(subs(C_sym,x,x0));

Dy = zeros(3,2);
%% Rumore sensori IMU

R = diag([
0.25^2;
0.1^2;
0.1^2
]);

Bw = [
0 0;
1/J_alpha_nom 0;
0 0;
0 1/Jbeta_nom
];


A = A_nom;

B = [B_nom Bw];

Cz = [
1 0 0 0;
0 0 1 0;
];

Dzu=zeros(2,2);

Dzw=zeros(2,2);

Dyu = zeros(3,2);

Dyw = zeros(size(Cy,1),size(Bw,2));

Cgen = [
Cz;
Cy
];

Dgen = [
    Dzu Dzw;
    Dyu Dyw
];

Pgen = ss(A,B,Cgen,Dgen);


%% LINEARIZZAZIONE INCERTA

% Creo direttamente matrici uncertain

A_unc = A_nom + 0*J_alpha_u;
B_unc = B_nom + 0*J_alpha_u;


Jbeta_u = Jy_u*sin(alpha_0)^2 + ...
          (Jz_u+m_u*l_u^2)*cos(alpha_0)^2 + ...
          p.I_b;


% Pitch
A_unc(2,1) = (-p.m*p.g*l_u*cos(alpha_0))/J_alpha_u;

A_unc(2,2) = -p.c_alpha/J_alpha_u;


B_unc(2,1) = l_u*cos(beta_0)/J_alpha_u;

B_unc(2,2) = eps_p_u*l_u*sin(beta_0)/J_alpha_u;



% Yaw
A_unc(4,4) = -p.c_beta/Jbeta_u;


B_unc(4,1) = eps_y_u*l_u*sin(alpha_0)/Jbeta_u;

B_unc(4,2) = l_u*cos(alpha_0)/Jbeta_u;
%% MODELLO NOMINALE

P_nom = ss(A_nom,B_nom,Cy,Dy);



%% MODELLO INCERTO

P_unc = uss(A_unc,B_unc,Cy,Dy);
P_unc.Uncertainty

disp('Modello nominale:')
P_nom

disp('impianto generalizzato:')
Pgen

disp('Modello incerto:')
P_unc

A_nom
B_nom

eig(A_nom)

tzero(P_nom)




rank(ctrb(A_nom,B_nom))
rank(obsv(A_nom,Cy))

% Inclusione della dinamica incerta degli attuatori (Utile per analisi successive)
P_full_nom = P_nom*blkdiag(G_act_nom,G_act_nom);

P_full_unc = P_unc*blkdiag(G_actuator_unc,G_actuator_unc);
P_full_unc.Uncertainty
disp('Matrici A e B calcolate. Impianto LTI incerto P_full creato con successo.');

