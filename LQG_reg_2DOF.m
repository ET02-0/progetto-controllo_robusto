
%% ==========================================
%  FASE 2: SINTESI LQG
%  Plant + Attuatori
%  (senza ritardo di Padé nella sintesi)
% ==========================================

disp('==============================================')
disp(' SINTESI LQG ESTESA')
disp('==============================================')



%% -------------------------------------------------
% 2) Modello attuatori
%    secondo ordine wn=30 rad/s zeta=0.7
% --------------------------------------------------

[A_act,B_act,C_act,D_act] = ssdata(G_act_2nd);

Actuator = ss(A_act,B_act,C_act,D_act);

% due attuatori indipendenti
Actuators_MIMO = blkdiag(Actuator,Actuator);


%% -------------------------------------------------
% 3) Connessione:
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

P_sintesi = P_nom*Actuators_MIMO;


[A_ext,B_ext,C_ext,D_ext]=ssdata(P_sintesi);

n_ext = size(A_ext,1);
m     = size(B_ext,2);
c     = size(C_ext,1);

M = [A_ext B_ext;
     -C_ext  zeros(c,m)];

disp(rank(M))
disp(size(M,1))
tzero(P_sintesi)




disp('Dimensioni modello esteso:')
fprintf('Stati   = %d\n',n_ext)
fprintf('Ingressi= %d\n',m)
fprintf('Uscite  = %d\n',c)


if n_ext ~= 8
    error('Errore: il modello esteso non ha 8 stati.')
end



%% -------------------------------------------------
% 4) Analisi modello esteso
% --------------------------------------------------

disp('Poli modello esteso:')
disp(eig(A_ext))


figure
pzmap(P_sintesi)
grid on
title('Poli modello esteso')

%% Controllabilità modello esteso

rank_ctrb = rank(ctrb(A_ext,B_ext));

fprintf('\nControllabilita LQG: %d/%d\n',...
    rank_ctrb,n_ext);


rank_obs = rank(obsv(A_ext,C_ext));

fprintf('Osservabilita LQG: %d/%d\n',...
    rank_obs,n_ext);


%% -------------------------------------------------
% 8) Pesi LQR (Regola di Bryson + Smorzamento)
% --------------------------------------------------
% Stato x = [alpha, alpha_dot, beta, beta_dot, x_act1, x_act1_dot, x_act2, x_act2_dot]

%% Esempio di tuning per smorzare le oscillazioni nel regolatore a 8 stati
%% Pesi LQR tuning smorzamento

%% Tuning maggiore smorzamento pitch

%% Tuning più smorzato

Q_heli = diag([800, 20, 10000, 500]);

Q_act = diag([5, 0.5, 5, 0.5]);

Q_lqr = blkdiag(Q_heli, Q_act);

R_lqr = diag([0.5 0.25]);
%% -------------------------------------------------
% 9) Calcolo LQR
% --------------------------------------------------

K_lqr = lqr(A_ext,B_ext,Q_lqr,R_lqr);


disp('Guadagni LQR')
disp(K_lqr)


Acl = A_ext-B_ext*K_lqr;

disp('Poli sistema chiuso LQR:')
disp(eig(Acl))



disp('Guadagno LQR calcolato')



%% -------------------------------------------------
% 10) Filtro di Kalman
% --------------------------------------------------


Gk = B_ext;

W = diag([
    0.001
    0.001
]);

V = sensor.R;

Qn = W;
Rn = V;

Estimator = ss(...
    A_ext,...
    Gk,...
    C_ext,...
    zeros(c,m));


[~,Ke,~] = kalman(Estimator,Qn,Rn);



disp('Filtro Kalman calcolato')



%% ==========================================
% Costruzione VERO controllore LQG a 2-DOF (Senza Integratore)
% ==========================================

% 1) Precompensatore Statico (Giusto, lo teniamo)
Acl = A_ext - B_ext*K_lqr;
C_track = [
1 0 0 0 0 0 0 0;
0 0 1 0 0 0 0 0
];
Kr = -inv(C_track*(Acl\B_ext));
disp('Precompensatore Kr:');
disp(Kr);

% 2) Matrici di stato del Controllore 2-DOF
% Lo stato è xc = x_hat
Ac_ctrl = A_ext - B_ext*K_lqr - Ke*C_ext;

% Il controllore ora ha 4 ingressi: i 2 riferimenti (r) e le 2 misure (y)
% L'ingresso r entra nell'osservatore moltiplicato per B_ext*Kr
% L'ingresso y entra nell'osservatore moltiplicato per Ke
Bc_ctrl = [B_ext*Kr, Ke]; 

% L'uscita del controllore è la u totale: u = -K_lqr*x_hat + Kr*r
Cc_ctrl = -K_lqr;
Dc_ctrl = [Kr, zeros(m,c)]; % Il termine Kr*r passa direttamente all'uscita

% Creazione del blocco
K_lqg_2dof = ss(Ac_ctrl, Bc_ctrl, Cc_ctrl, Dc_ctrl);
save('LQG_Controllers.mat','K_lqg_2dof','-append');
K_lqg_2dof.InputName = {
'r_alpha',...
'r_beta',...
'y_acc',...
'm_x',...
'm_y'};
K_lqg_2dof.OutputName = {'u1', 'u2'};

disp('Controllore LQG 2-DOF creato con successo');

%% ==========================================
% Modello di monitoraggio alpha beta
% ==========================================

%% Plant per monitoraggio alpha beta + feedback LQG

C_monitor_track = [
    C_ext;
    1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0
];

D_monitor_track = zeros(5,2);

P_monitor = ss(A_ext,B_ext,C_monitor_track,D_monitor_track);


P_monitor.InputName = {'u1','u2'};

P_monitor.OutputName = {
    'y_acc',...
    'm_x',...
    'm_y',...
    'alpha',...
    'beta'};
CL_monitor = connect(P_monitor,...
                     K_lqg_2dof,...
                     {'r_alpha','r_beta'},...
                     {'alpha','beta'});
CL_monitor.OutputName

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
    zeros(4,2)
];

B_dist = [B_ext Bd_ext];

% Uscite:
% 1-3 -> sensori
% 4   -> delta_alpha
% 5   -> delta_beta

C_dist = [
    C_ext;
    1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0
];

D_dist = zeros(5,4);

P_dist = ss( ...
    A_ext, ...
    B_dist, ...
    C_dist, ...
    D_dist);

P_dist.InputName = {
    'u1'
    'u2'
    'd_alpha'
    'd_beta'
};

P_dist.OutputName = {
    'y_acc'
    'm_x'
    'm_y'
    'delta_alpha'
    'delta_beta'
};

% Collegamento closed-loop:
% ingressi esterni = r_alpha, r_beta, d_alpha, d_beta
%
% uscite richieste = delta_alpha, delta_beta

CL_dist = connect( ...
    P_dist, ...
    K_lqg_2dof, ...
    {'r_alpha','r_beta','d_alpha','d_beta'}, ...
    {'delta_alpha','delta_beta'});

%% Trasferimento d_beta -> delta_beta

G_dBeta_beta = CL_dist(2,4);

disp('==============================================')
disp(' TRASFERIMENTO d_beta -> delta_beta')
disp('==============================================')

disp(G_dBeta_beta)

dcgain_dBeta_beta = dcgain(G_dBeta_beta);

fprintf('\nGuadagno statico d_beta -> delta_beta = %.8f\n', ...
    dcgain_dBeta_beta);

%% Risposta a gradino del disturbo beta

figure

step(G_dBeta_beta)

grid on

xlabel('Tempo [s]')
ylabel('\Delta\beta [rad]')
title('Reiezione disturbo: d_\beta \rightarrow \Delta\beta')

%% Risposta al disturbo con ampiezza usata in Simulink

d_beta_test = 2e-3;

figure

step(d_beta_test * G_dBeta_beta)

grid on

xlabel('Tempo [s]')
ylabel('\Delta\beta [rad]')
title('Risposta a d_\beta = 2e-3 N m')

fprintf('\nRisposta statica prevista per d_beta = %.4g N m:\n', ...
    d_beta_test);

fprintf('delta_beta_ss = %.8f rad (%.4f deg)\n', ...
    dcgain_dBeta_beta*d_beta_test, ...
    rad2deg(dcgain_dBeta_beta*d_beta_test));

V = sensor.R;

Vinv = diag(1./diag(V));

C_angle_meas = C_ext(:,[1 3]);

HalphaBeta = ...
    (C_angle_meas.'*Vinv*C_angle_meas) \ ...
    (C_angle_meas.'*Vinv);
%% ==========================================
% Connessione LQG regolatore
% ==========================================


P_full_nom.InputName = {'u1','u2'};
P_full_nom.OutputName = {
'y_acc',...
'm_x',...
'm_y'};


CL_nom = connect(P_full_nom,...
                 K_lqg_2dof,...
                 {'r_alpha','r_beta'},...
                 {'y_acc','m_x','m_y'});




disp('==============================================')
disp(' ANALISI STABILITA')
disp('==============================================')

% Poli sistema chiuso nominale
poli_CL = pole(CL_nom);

disp('Poli closed-loop:')
disp(poli_CL)


if all(real(poli_CL)<0)
    disp('Sistema nominale stabile')
else
    disp('ATTENZIONE: sistema nominale INSTABILE')
end


% Margine di stabilità classico
figure
pzmap(CL_nom)
grid on
title('Poli Closed Loop LQG')


% Risposta al gradino
figure
step(CL_nom)

grid on
title('Risposta LQG nominale')

P_full_unc.InputName = {'u1','u2'};
P_full_unc.OutputName = {
'y_acc',...
'm_x',...
'm_y'};

K_lqg_2dof.InputName = {
'r_alpha',...
'r_beta',...
'y_acc',...
'm_x',...
'm_y'};
K_lqg_2dof.OutputName = {'u1','u2'};


CL_unc = connect(P_full_unc,...
                 K_lqg_2dof,...
                 {'r_alpha','r_beta'},...
                 {'y_acc','m_x','m_y'});
%% ==========================================
% Analisi robustezza
% ==========================================

disp('==============================================')
disp(' ANALISI ROBUSTEZZA')
disp('==============================================')


% Robust stability margin
stab_margin = robuststab(CL_unc);

disp('Margine robustezza:')
disp(stab_margin)


% Diagramma valori singolari
omega = logspace(-2,3,500);

figure
sigma(CL_unc,omega)
grid on
title('Singular Values Closed Loop incerto')

%% Robust performance
%% Analisi prestazioni

[y_angle,t_angle] = step(CL_monitor);

figure

subplot(2,2,1)
plot(t_angle,y_angle(:,1,1))
grid on
title('\alpha / r_\alpha')

subplot(2,2,2)
plot(t_angle,y_angle(:,1,2))
grid on
title('\alpha / r_\beta')

subplot(2,2,3)
plot(t_angle,y_angle(:,2,1))
grid on
title('\beta / r_\alpha')

subplot(2,2,4)
plot(t_angle,y_angle(:,2,2))
grid on
title('\beta / r_\beta')

%{
%% ==========================================
% TEST DI PROGETTO
% r_alpha = 0.1 rad
% r_beta  = 0 rad
% ==========================================

t = linspace(0,20,2001);

r = zeros(length(t),2);
r(:,1) = 0.1;
r(:,2) = 0;

[y_lsim,t_lsim] = lsim(CL_monitor,r,t);

alpha_resp = y_lsim(:,1);
beta_resp  = y_lsim(:,2);

figure
plot(t_lsim,alpha_resp,'LineWidth',1.5)
hold on
plot(t_lsim,beta_resp,'LineWidth',1.5)
grid on
xlabel('Tempo [s]')
ylabel('Angolo [rad]')
legend('\alpha','\beta','Location','best')
title('LQG 2DOF senza integratore - r_\alpha=0.1 rad, r_\beta=0')

%% Metriche del caso di progetto

info_alpha = stepinfo(alpha_resp,t_lsim);

fprintf('\n==============================================\n')
fprintf(' PRESTAZIONI LQG SENZA INTEGRATORE\n')
fprintf('==============================================\n')

fprintf('Overshoot alpha      = %.4f %%\n',...
    info_alpha.Overshoot);

fprintf('Settling time alpha  = %.4f s\n',...
    info_alpha.SettlingTime);

fprintf('Rise time alpha      = %.4f s\n',...
    info_alpha.RiseTime);

fprintf('Errore finale alpha  = %.6f rad\n',...
    abs(alpha_resp(end)-0.1));

fprintf('Picco |alpha|        = %.6f rad\n',...
    max(abs(alpha_resp)));

fprintf('Picco |beta|         = %.6f rad\n',...
    max(abs(beta_resp)));

fprintf('RMS beta             = %.6f rad\n',...
    rms(beta_resp));

fprintf('Errore finale beta   = %.6f rad\n',...
    abs(beta_resp(end)));

%}

gain_aa = dcgain(CL_monitor(1,1));
gain_bb = dcgain(CL_monitor(2,2));
gain_ab = dcgain(CL_monitor(1,2));
gain_ba = dcgain(CL_monitor(2,1));

fprintf('\nGuadagni statici:\n')
fprintf('Gaa = %.6f\n',gain_aa)
fprintf('Gbb = %.6f\n',gain_bb)
fprintf('Gab = %.6e\n',gain_ab)
fprintf('Gba = %.6e\n',gain_ba)

fprintf('Errore statico alpha = %.6f %%\n',100*abs(1-gain_aa))
fprintf('Errore statico beta  = %.6f %%\n',100*abs(1-gain_bb))
%{
% 1. Creo l'oggetto opzioni per l'analisi Worst-Case e forzo l'algoritmo 'a' (Advanced)
opt = wcOptions('MussvOptions', 'a');

% 2. Calcolo il Worst-Case Gain passando le opzioni
[wcg, wcu] = wcgain(CL_unc, opt);

disp('Worst-case gain (Algoritmo Avanzato):')
disp(wcg)

robustperf(CL_unc)
%}