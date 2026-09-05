% =========================================================================
% Script di Inizializzazione e Linearizzazione: Elicottero 2DoF
% Progetto di Controllo Robusto (Sintesi Robusta su 4 stati)
% =========================================================================
clear all; close all; clc;
disp('Configurazione parametri Elicottero 2DoF (Controllo Robusto)...');
which uss -all

%% 1. PARAMETRI NOMINALI (Perfettamente noti)
p.J_y = 0.00023; 
p.J_z = 0.00364; 
p.I_b = 0.00023; 
p.m = 0.2; 
p.l_nom = 0.2; 
p.c_alpha = 0.01; 
p.c_beta = 0.01;
p.g = 9.81;

%% 2. PUNTO DI EQUILIBRIO (Trim Point)

alpha_0 = deg2rad(10);
beta_0  = 0;

% --- Condizioni iniziali per i blocchi Integratore di Simulink ---
q0    = [alpha_0; beta_0]; 
qdot0 = [0; 0];

%% PARAMETRI NOMINALI PER LQG

J_alpha_nom = 0.012;
l_nom       = 0.2;
eps_p_nom   = 0.1;
eps_y_nom   = 0.1;
tau_d_nom   = 0.02;
Jbeta_nom = p.J_y*sin(alpha_0)^2 + ...
            (p.J_z+p.m*l_nom^2)*cos(alpha_0)^2 + ...
            p.I_b;

% Allineamento nomi per la MATLAB Function di Simulink
p.l = p.l_nom;
p.eps_p = eps_p_nom;
p.eps_y = eps_y_nom;
p.J_alpha = J_alpha_nom;


%% EQUILIBRIO NOMINALE

E0_nom = [
    cos(beta_0),              eps_p_nom*sin(beta_0);
    eps_y_nom*sin(alpha_0),   cos(alpha_0)
];

u0_nom = E0_nom \ [
    p.m*p.g*sin(alpha_0);
    0
];

F1_0 = u0_nom(1);
F2_0 = u0_nom(2);
% CALCOLO DELLE MISURE DI EQUILIBRIO (TRIM DEI SENSORI)
y_acc_eq = -p.g * sin(alpha_0);
mx_eq    = cos(alpha_0) * cos(beta_0);
my_eq    = -sin(beta_0);

y0 = [y_acc_eq; mx_eq; my_eq]; % Questo è il tuo y_eq



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

%% ========================================================================
% 4. STRUTTURE PER SIMULINK (Attuatori, Sensori, Disturbi, Riferimenti)
% ========================================================================

% --- ATTUATORI (act) ---
% Raggruppiamo i parametri dinamici scelti in precedenza (omega_n, zeta, tau_d_nom)
act.wn1 = omega_n;       
act.wn2 = omega_n;       
act.zeta1 = zeta;        
act.zeta2 = zeta;        
act.td1 = tau_d_nom;     
act.td2 = tau_d_nom;     

% Saturazioni fisiche per la validazione della RAS
% (Centrate in modo da consentire variazioni attorno a F1_0 e F2_0)
act.deltaF1_min = -0.5;  
act.deltaF1_max =  2.5;  
act.deltaF2_min = -2.5;  
act.deltaF2_max =  2.5;  

% Rumore di attuazione
act.Ts_noise = 0.01;    % Campionamento del rumore attuatori
act.sigma_F1 = 0.01;     % [N]
act.sigma_F2 = 0.01;     % [N]
act.noisePower_F1 = act.sigma_F1^2 * act.Ts_noise;
act.noisePower_F2 = act.sigma_F2^2 * act.Ts_noise;
act.noiseEnable = 1;     % Mettere a 0 per disattivare il rumore attuatori

% --- SENSORI IMU MPU-9250 (sensor) ---
sensor.seed = 12345;     
sensor.noiseEnable = 1;  % Mettere a 0 per la valutazione teorica della RAS

% Accelerometro (Campionato ad alta frequenza)
sensor.acc.fs = 1000;              
sensor.acc.Ts = 1/sensor.acc.fs;   % dt = 0.001 s
sensor.acc.sigma = 0.25;           
sensor.acc.var = sensor.acc.sigma^2; 

% Magnetometro (Campionato a bassa frequenza)
sensor.mag.fs = 100;               
sensor.mag.Ts = 1/sensor.mag.fs;   % T_cmag = 0.01 s
sensor.mag.sigma = 0.1;         
sensor.mag.var = sensor.mag.sigma^2; 

sensor.acc.noisePower = sensor.acc.sigma^2 * sensor.acc.Ts;
sensor.mag.noisePower = sensor.mag.sigma^2 * sensor.mag.Ts;

sensor.R = diag([
    sensor.acc.sigma^2
    sensor.mag.sigma^2
    sensor.mag.sigma^2
]);


% --- DISTURBI AERODINAMICI (aero) ---
% Come richiesto dalla Traccia 3 per la verifica della reiezione ai disturbi
aero.enable = 1;                  
aero.alpha.time      = 10;         % Istante ingresso disturbo su pitch [s]
aero.alpha.amplitude = 5e-3;      % Ampiezza disturbo [N*m]
aero.beta.time       = 10;        % Istante ingresso disturbo su yaw [s]
aero.beta.amplitude  = 2e-3;      % Ampiezza disturbo [N*m]

% --- RIFERIMENTI (ref) ---
% Impostiamo un test di gradino su alpha mantenendo beta all'equilibrio nominale
ref.alpha.initial = 0;
ref.alpha.final   = deg2rad(3);
ref.alpha.time    = 1;
ref.beta.initial  = 0;
ref.beta.final    = 0;
ref.beta.time     = 1;

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

umax = 5;

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

%% EQUILIBRIO INCERTO

E0_unc = [
    cos(beta_0),              eps_p_u*sin(beta_0);
    eps_y_u*sin(alpha_0),     cos(alpha_0)
];

u0_unc = E0_unc \ [
    m_u*p.g*sin(alpha_0);
    0
];

F1_0_unc = u0_unc(1);
F2_0_unc = u0_unc(2);

%% LINEARIZZAZIONE INCERTA

Jbeta_u = ...
    Jy_u*sin(alpha_0)^2 + ...
    (Jz_u + m_u*l_u^2)*cos(alpha_0)^2 + ...
    p.I_b;

% ---------- MATRICE A INCERTA ----------

a21_u = ...
    -(m_u*p.g*l_u*cos(alpha_0))/J_alpha_u;

a22_u = ...
    -p.c_alpha/J_alpha_u;

a23_u = ...
    l_u*( ...
    -F1_0_unc*sin(beta_0) + ...
    eps_p_u*F2_0_unc*cos(beta_0)) / J_alpha_u;

a41_u = ...
    l_u*( ...
    -F2_0_unc*sin(alpha_0) + ...
    eps_y_u*F1_0_unc*cos(alpha_0)) / Jbeta_u;

a44_u = ...
    -p.c_beta/Jbeta_u;

A_unc = [
    0      1      0      0;
    a21_u a22_u  a23_u  0;
    0      0      0      1;
    a41_u 0      0      a44_u
];

% ---------- MATRICE B INCERTA ----------

b21_u = l_u*cos(beta_0)/J_alpha_u;

b22_u = eps_p_u*l_u*sin(beta_0)/J_alpha_u;

b41_u = eps_y_u*l_u*sin(alpha_0)/Jbeta_u;

b42_u = l_u*cos(alpha_0)/Jbeta_u;

B_unc = [
    0      0;
    b21_u b22_u;
    0      0;
    b41_u b42_u
];
%% MODELLO NOMINALE
P_nom = ss(A_nom,B_nom,Cy,Dy);

P_nom.StateName = {
    'delta_alpha'
    'delta_alpha_dot'
    'delta_beta'
    'delta_beta_dot'
};

P_nom.InputName = {
    'delta_F1'
    'delta_F2'
};

P_nom.OutputName = {
    'delta_y_acc'
    'delta_mx'
    'delta_my'
};

%% MODELLO NOMINALE CON DISTURBI AERODINAMICI

Bd_nom = [
    0                 0;
    1/J_alpha_nom     0;
    0                 0;
    0                 1/Jbeta_nom
];

B_ext_nom = [B_nom Bd_nom];

P_nom_ext = ss( ...
    A_nom, ...
    B_ext_nom, ...
    Cy, ...
    zeros(3,4));

P_nom_ext.StateName = P_nom.StateName;

P_nom_ext.InputName = {
    'delta_F1'
    'delta_F2'
    'd_alpha'
    'd_beta'
};

P_nom_ext.OutputName = P_nom.OutputName;

%% Attuatori
%% -------------------------------------------------
% Connessione:
%
%       u
%       |
%       v
%   Attuatori
%       |
%       v
%   Elicottero
%       |
%       v
%       y
%
% --------------------------------------------------


[A_act,B_act,C_act,D_act] = ssdata(G_act_2nd);

Actuator = ss(A_act,B_act,C_act,D_act);

Actuators_MIMO = blkdiag(Actuator,Actuator);


%% Plant esteso

P_ext = P_nom*Actuators_MIMO;

[A_ext,B_ext,C_ext,D_ext]=ssdata(P_ext);

M = [A_ext B_ext;
     -C_ext zeros(size(C_ext,1),size(B_ext,2))];

V = sensor.R;

Vinv = diag(1./diag(V));

C_angle_meas = C_ext(:,[1 3]);

HalphaBeta = ...
    (C_angle_meas.'*Vinv*C_angle_meas) \ ...
    (C_angle_meas.'*Vinv);

%% MODELLO INCERTO
P_unc = uss(A_unc,B_unc,Cy,Dy);

P_unc.StateName = P_nom.StateName;
P_unc.InputName = P_nom.InputName;
P_unc.OutputName = P_nom.OutputName;


%% DISTURBI AERODINAMICI
Bd_unc = [
    0             0;
    1/J_alpha_u  0;
    0             0;
    0             1/Jbeta_u
];

B_ext_unc = [B_unc Bd_unc];

P_unc_ext = uss( ...
    A_unc, ...
    B_ext_unc, ...
    Cy, ...
    zeros(3,4));

P_unc_ext.StateName = P_unc.StateName;

P_unc_ext.InputName = {
    'delta_F1'
    'delta_F2'
    'd_alpha'
    'd_beta'
};

P_unc_ext.OutputName = {
    'delta_y_acc'
    'delta_mx'
    'delta_my'
};


%% MODELLO INCERTO CON USCITE ANGOLARI
C_angles = [
    1 0 0 0;
    0 0 1 0
];

P_unc_angles = uss( ...
    A_unc, ...
    B_unc, ...
    C_angles, ...
    zeros(2,2));

P_unc_angles.StateName = P_unc.StateName;

P_unc_angles.InputName = P_unc.InputName;

P_unc_angles.OutputName = {
    'delta_alpha'
    'delta_beta'
};


%% ATTUATORI
P_full_nom = P_nom*blkdiag(G_act_2nd,G_act_2nd);
P_full_unc = P_unc*blkdiag(G_actuator_unc,G_actuator_unc);


P_full_unc.Uncertainty

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

disp('Matrici A e B calcolate. Impianto LTI incerto P_full creato con successo.');

