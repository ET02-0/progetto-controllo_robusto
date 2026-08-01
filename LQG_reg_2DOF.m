
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

M = [A_ext B_ext;
     -C_ext  zeros(p,m)];

disp(rank(M))
disp(size(M,1))
tzero(P_sintesi)

n_ext = size(A_ext,1);
m     = size(B_ext,2);
p     = size(C_ext,1);


disp('Dimensioni modello esteso:')
fprintf('Stati   = %d\n',n_ext)
fprintf('Ingressi= %d\n',m)
fprintf('Uscite  = %d\n',p)


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

Q_heli = diag([150,350,150,250]);

Q_act = diag([5,0.5,5,0.5]);

Q_lqr = blkdiag(Q_heli,Q_act);

R_lqr = 3*eye(2);
%% -------------------------------------------------
% 9) Calcolo LQR
% --------------------------------------------------

K_lqr = lqr(A_ext,B_ext,Q_lqr,R_lqr);
umax = 5; % esempio Nm

disp('Guadagni LQR')
disp(K_lqr)


Acl = A_ext-B_ext*K_lqr;

disp('Poli sistema chiuso LQR:')
disp(eig(Acl))



disp('Guadagno LQR calcolato')



%% -------------------------------------------------
% 10) Filtro di Kalman
% --------------------------------------------------

% Assumiamo che il rumore di processo entri come disturbo di attuazione (vento)
% Rumore di processo modellato come disturbo sugli ingressi
% (variazioni di spinta, dinamica non modellata degli attuatori)
Gk = B_ext;

Qn = 5e-3*eye(m);
Rn = R;

Estimator = ss(...
    A_ext,...
    Gk,...
    C_ext,...
    zeros(p,m));


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
Dc_ctrl = [Kr, zeros(m,p)]; % Il termine Kr*r passa direttamente all'uscita

% Creazione del blocco
K_lqg_2dof = ss(Ac_ctrl, Bc_ctrl, Cc_ctrl, Dc_ctrl);
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

%% Metriche tracking diagonali

info_aa = stepinfo(y_angle(:,1,1),t_angle);
info_bb = stepinfo(y_angle(:,2,2),t_angle);

disp('Overshoot alpha tracking:')
disp(info_aa.Overshoot)

disp('Overshoot beta tracking:')
disp(info_bb.Overshoot)

disp('Guadagno statico alpha-alpha')
disp(y_angle(end,1,1))

disp('Guadagno statico beta-beta')
disp(y_angle(end,2,2))

disp('Accoppiamento alpha-beta')
disp(y_angle(end,1,2))

disp('Accoppiamento beta-alpha')
disp(y_angle(end,2,1))

disp('Settling time alpha:')
disp(info_aa.SettlingTime)

disp('Settling time beta:')
disp(info_bb.SettlingTime)


disp('Rise time alpha:')
disp(info_aa.RiseTime)

disp('Rise time beta:')
disp(info_bb.RiseTime)

% 1. Creo l'oggetto opzioni per l'analisi Worst-Case e forzo l'algoritmo 'a' (Advanced)
opt = wcOptions('MussvOptions', 'a');

% 2. Calcolo il Worst-Case Gain passando le opzioni
[wcg, wcu] = wcgain(CL_unc, opt);

disp('Worst-case gain (Algoritmo Avanzato):')
disp(wcg)

robustperf(CL_unc)