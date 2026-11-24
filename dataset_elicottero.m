% =========================================================================
% Script di Inizializzazione e Linearizzazione: Elicottero 2DoF
% Progetto di Controllo Robusto (Sintesi Robusta su 4 stati)
% =========================================================================
clear; close all; clc;
disp('Configurazione parametri Elicottero 2DoF (Controllo Robusto)...');
which uss -all

% ================================================================
% IMPOSTAZIONE GRAFICA PER LA RELAZIONE
% ================================================================

set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesColor','w', ...
    'defaultAxesXColor','k', ...
    'defaultAxesYColor','k', ...
    'defaultAxesZColor','k', ...
    'defaultAxesGridColor',[0.75 0.75 0.75], ...
    'defaultTextColor','k', ...
    'defaultAxesBox','on', ...
    'defaultAxesFontSize',11);

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
beta_0  = deg2rad(15);

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

y0 = [y_acc_eq; mx_eq; my_eq]; % y_eq



%% PARAMETRI INCERTI PER ROBUSTEZZA

J_alpha_u = ureal('J_alpha',0.012,'Percentage',10);
l_u       = ureal('l',0.2,'Percentage',5);
eps_p_u   = ureal('eps_p',0.1,'Percentage',20);
eps_y_u   = ureal('eps_y',0.1,'Percentage',20);
Jy_u = ureal('Jy',0.00023,'Percentage',5);
Jz_u = ureal('Jz',0.00364,'Percentage',5);
m_u  = ureal('m',0.2,'Percentage',5);
omega_n_u1 = ureal('omega_n1',30,'Percentage',10);
omega_n_u2 = ureal('omega_n2',30,'Percentage',10);
tau_d_u1   = ureal('tau_d1',0.02,'Percentage',20);
tau_d_u2   = ureal('tau_d2',0.02,'Percentage',20);

%% ============================================================
%  CASO SEPARATO:
%  omega_n e tau_d parametrici + residuo complesso Padé
% =============================================================

s = tf('s');

% Intervallo ammesso di tau_d
tau_vec = linspace(0.80*tau_d_nom, 1.20*tau_d_nom, 101);

% Griglia di frequenza per il fit di W_Pade
w = logspace(-2,4,2000);

% Famiglia del residuo:
%
% Residuo additivo:
% E_pade(s,tau) = exp(-tau*s) - G_pade(s,tau)


%% ============================================================
% ERRORE ADDITIVO DELLA APPROSSIMAZIONE PADE
% ============================================================

E_pade_array = [];

for k = 1:length(tau_vec)

    tau = tau_vec(k);

    % Ritardo esatto
    G_delay_exact = tf(1,1,'InputDelay',tau);

    % Approssimazione di Pade del primo ordine
    G_delay_pade = ...
        (1 - tau*s/2) / ...
        (1 + tau*s/2);

    % Residuo ADDITIVO: ritardo esatto meno Pade
    E_pade = G_delay_exact - G_delay_pade;

    if k == 1
        E_pade_array = E_pade;
    else
        E_pade_array = cat(3,E_pade_array,E_pade);
    end

end

% Conversione in FRD
E_pade_frd = frd(E_pade_array,w);

%% ============================================================
% INVILUPPO DELL'ERRORE ADDITIVO
% ============================================================

[resp_E,freq_R] = frdata(E_pade_frd);

resp_E = squeeze(resp_E);

% Modulo dell'errore additivo
err_pade = abs(resp_E);

% Inviluppo massimo sui valori di tau_d
err_env = max(err_pade,[],2);

% Evita lo zero esatto durante il fit
err_floor = 1e-8;
err_env_fit = max(err_env,err_floor);

%% ============================================================
% FIT DEL PESO W_Pade
% ============================================================

E_pade_env = frd( ...
    reshape(err_env_fit,1,1,[]), ...
    freq_R);

% Manteniamo il peso di ordine 2
order_pade = 4;

% Margine di copertura del 10%
Cfit.LowerBound = 1.10 * E_pade_env;
Cfit.UpperBound = Inf;

W_Pade = fitmagfrd( ...
    E_pade_env, ...
    order_pade, ...
    [], ...
    [], ...
    Cfit);


%% ================================================================
% VERIFICA W_Pade
% ================================================================

disp('============================================================');
disp(' VERIFICA W_Pade');
disp('============================================================');

W_Pade

disp('Poli W_Pade:');
pole(W_Pade)

mag_WPade = squeeze(abs(freqresp(W_Pade,freq_R)));

fprintf('NaN presenti in W_Pade: %d\n', ...
    any(isnan(mag_WPade)));

fprintf('Inf presenti in W_Pade: %d\n', ...
    any(isinf(mag_WPade)));

fprintf('Modulo minimo W_Pade: %.6e\n', ...
    min(mag_WPade));

fprintf('Modulo massimo W_Pade: %.6e\n', ...
    max(mag_WPade));

%% ================================================================
% VERIFICA NUMERICA DEL COVER ADDITIVO
% ================================================================

ratio_Pade = err_env ./ mag_WPade;

% Rapporto che include il margine richiesto del 10%
ratio_Pade_10 = 1.10 * err_env ./ mag_WPade;

[max_ratio,idx_max] = max(ratio_Pade);
max_ratio_10 = max(ratio_Pade_10);

fprintf('\n============================================================\n');
fprintf(' RISULTATO COVER ADDITIVO PADE\n');
fprintf('============================================================\n');

fprintf('Massimo errore additivo = %.6e\n', max(err_env));
fprintf('Massimo modulo W_Pade   = %.6e\n', max(mag_WPade));

fprintf('Massimo rapporto errore/W = %.6f\n', max_ratio);
fprintf('Massimo rapporto 1.10*errore/W = %.6f\n', max_ratio_10);

fprintf('Frequenza del massimo rapporto = %.6f rad/s\n', ...
    freq_R(idx_max));

fprintf('Punti fuori copertura base = %d / %d\n', ...
    sum(ratio_Pade > 1 + 1e-6), numel(ratio_Pade));

fprintf('Punti che violano il margine del 10%% = %d / %d\n', ...
    sum(ratio_Pade_10 > 1 + 1e-6), numel(ratio_Pade_10));
W_hf = abs(squeeze(freqresp(W_Pade,1e6)));

fprintf('Modulo W_Pade a 1e6 rad/s = %.6f\n', W_hf);

fprintf('Rapporto massimo errore/peso: %.6f\n', ...
    max(err_env ./ mag_WPade));

fprintf('Rapporto massimo con margine 10%%: %.6f\n', ...
    max(1.10 * err_env ./ mag_WPade));

%% ================================================================
% GRAFICO DEL COVER
% ================================================================

figure('Name','Cover errore Padé');

semilogx( ...
    freq_R, ...
    err_env, ...
    'LineWidth',1.5);

hold on;

semilogx( ...
    freq_R, ...
    mag_WPade, ...
    'LineWidth',1.5);

grid on;

xlabel('\omega [rad/s]');
ylabel('Modulo');

title('Inviluppo errore Padé e W_{Pade}');

legend( ...
    '|G_{ritardo}-G_{Pade}| massimo', ...
    '|W_{Pade}|', ...
    'Location','best');

%% ================================================================
% INCERTEZZA COMPLESSA DINAMICA
% ================================================================

Delta_Pade1 = ultidyn('Delta_Pade1',[1 1]);
Delta_Pade2 = ultidyn('Delta_Pade2',[1 1]);

%% 3. DINAMICA DEGLI ATTUATORI

omega_n = 30;
zeta    = 0.7;

% ================================================================
% ATTUATORE NOMINALE
% ================================================================

num_2nd = omega_n^2;
den_2nd = [1, 2*zeta*omega_n, omega_n^2];

G_act_2nd = tf(num_2nd,den_2nd);


% Ritardo nominale tramite Padé del primo ordine

num_pade_nom = [-tau_d_nom/2 1];
den_pade_nom = [ tau_d_nom/2 1];

G_delay_pade_nom = tf(num_pade_nom,den_pade_nom);

G_act_nom = G_act_2nd * G_delay_pade_nom;


% ================================================================
% ATTUATORE INCERTO
%
% omega_n : +/-10%
% tau_d   : +/-20%
% zeta    : nominale
% ================================================================

s = tf('s');

% Attuatore 1
G_act_2nd_unc1 = omega_n_u1^2/(s^2 + 2*zeta*omega_n_u1*s + omega_n_u1^2);
G_delay_pade_unc1 = tf([-tau_d_u1/2 1],[tau_d_u1/2 1]);
G_actuator_unc1 = G_act_2nd_unc1*(G_delay_pade_unc1 + W_Pade*Delta_Pade1);

% Attuatore 2
G_act_2nd_unc2 = omega_n_u2^2/(s^2 + 2*zeta*omega_n_u2*s + omega_n_u2^2);
G_delay_pade_unc2 = tf([-tau_d_u2/2 1],[tau_d_u2/2 1]);
G_actuator_unc2 = G_act_2nd_unc2*(G_delay_pade_unc2 + W_Pade*Delta_Pade2);

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
act.deltaF1_min = -2.5;  
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

%% ==========================================
% ANALISI REIEZIONE DISTURBO
% Trasferimento d_beta -> delta_beta
% LQG 2-DOF senza integratore
% ==========================================

% Gli attuatori sono già inclusi in A_ext/B_ext.
% Aggiungiamo al modello esteso gli ingressi di disturbo
% d_alpha e d_beta direttamente sulla dinamica dell'elicottero.

Bd_ext = [
    Bd_nom
    zeros(6,2)
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

C_tracking = [
    1 0 0 0;
    0 0 1 0
    ];

C_out = [
    P_nom_ext.C;
    C_tracking
    ];

D_out = [
    P_nom_ext.D;
    zeros(2,size(P_nom_ext.B,2))
    ];


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


Actuator = ss(G_act_nom);

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

P_full_unc = P_unc*blkdiag(G_actuator_unc1, G_actuator_unc2);

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

