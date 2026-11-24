% =========================================================================
% Script di Inizializzazione e Linearizzazione: Elicottero 2DoF
% Progetto di Controllo Robusto
% LQG 2-DOF con integratore - modello esteso a 10 stati
% =========================================================================

%close all
clc

disp('SINTESI LQG 2DOF CON INTEGRATORE')


disp(rank(M))
disp(size(M,1))
tzero(P_ext)

n_ext = size(A_ext,1);
m     = size(B_ext,2);
c     = size(C_ext,1);


disp('Dimensioni modello esteso:')
fprintf('Stati   = %d\n',n_ext)
fprintf('Ingressi= %d\n',m)
fprintf('Uscite  = %d\n',c)


if n_ext ~= 10
    error('Errore: il modello esteso non ha 10 stati.')
end



%% -------------------------------------------------
% 4) Analisi modello esteso
% --------------------------------------------------

disp('Poli modello esteso:')
disp(eig(A_ext))


figure
pzmap(P_ext)
grid on
title('Poli e zeri modello esteso')


%% Integratori su pitch (alpha) e yaw (beta)
C_int = [
    1 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0
];

% Verifica coerenza tra uscita misurata e uscita da integrare
err_H = norm(HalphaBeta*C_ext - C_int, 'fro');

fprintf('Errore HalphaBeta*C_ext - C_int = %.3e\n', err_H);

ni = size(C_int,1);         % Ora ni = 2

Aa = [
A_ext zeros(n_ext,ni);
-C_int zeros(ni,ni)
];


Ba = [
B_ext;
zeros(ni,m)
];



%% -------------------------------------------------
% 6) Verifica controllabilità sistema aumentato - PBH
% -------------------------------------------------

lambda_a = eig(Aa);

controllabile_a = true;

for k = 1:length(lambda_a)

    Mcont_a = [
        lambda_a(k)*eye(size(Aa)) - Aa, ...
        Ba
    ];

    r = rank(Mcont_a);

    fprintf('PBH controllabilità, polo %d: rango = %d/%d\n', ...
        k, r, size(Aa,1));

    if r < size(Aa,1)
        controllabile_a = false;
    end
end

if controllabile_a
    disp('Il sistema aumentato è COMPLETAMENTE CONTROLLABILE (PBH).');
else
    warning('Il sistema aumentato NON è completamente controllabile (PBH).');
end


%% -------------------------------------------------
% 7) Verifica osservabilità - PBH
% -------------------------------------------------

lambda_obs = eig(A_ext);

osservabile = true;

for k = 1:length(lambda_obs)

    Mobs = [
        lambda_obs(k)*eye(n_ext) - A_ext;
        C_ext
    ];

    if rank(Mobs) < n_ext
        osservabile = false;
    end

    fprintf('PBH osservabilità, polo %d: rango = %d/%d\n', ...
        k, rank(Mobs), n_ext);
end

fprintf('\n');

if osservabile
    disp('Il modello esteso è COMPLETAMENTE OSSERVABILE (PBH).');
else
    warning('Il modello esteso NON è completamente osservabile (PBH).');
end

%% -------------------------------------------------
% 8) Pesi LQR
%
% Stato:
% x = [elicottero(4); attuatori(6)]
%
% [alpha, alphadot, beta, betadot,
%  act1_1, act1_2, act1_pade,
%  act2_1, act2_2, act2_pade]
%
% --------------------------------------------------
% =================================================================
% 1. Pesi LQR
% =================================================================


xmax = [deg2rad(10);   % alpha: 10°
        0.75;           % alpha_dot: 0.75 rad/s
        deg2rad(7);   % beta: 7°
        0.5  ];          % beta_dot: 0.5 rad/s

Q_heli_base = diag(1 ./ xmax.^2);

% Candidato B selezionato dalla ricerca
sa = 1;
sb = 1;

Q_heli = diag(diag(Q_heli_base) .* ...
              [sa; sa; sb; sb]);

deltaF1_max = 0.75;
deltaF2_max = 0.75;
% Matrice di uscita degli attuatori
[A_act, ~, C_act, D_act] = ssdata(Actuators_MIMO);

W_deltaF = diag(1 ./ [deltaF1_max; deltaF2_max].^2);
Q_act = C_act' * W_deltaF * C_act;


rhoR =   1.7783 ;  % Cambia questo valore per ogni test
R_lqr = rhoR * eye(2);



%% -------------------------------------------------
% 10) Filtro di Kalman
% --------------------------------------------------
Gk = [
    B_nom;
    zeros(6,2)
    ];

W = diag([
    30*act.sigma_F1^2
    15*act.sigma_F2^2
]);


V = sensor.R;

Qn = W;
Rn = V;



[Ke,~,~] = lqe(A_ext,Gk,C_ext,W,V);

disp('Filtro Kalman calcolato')



%% Selezione finale mediante tau_alpha e tau_beta


q_int_target = [140; 25];

% Valori massimi associati agli integratori di alpha e beta
xmax_int = xmax([1; 3]);

% Calcolo automatico dei due tau
tau_i = 1 ./ (xmax_int .* sqrt(q_int_target));

tau_alpha = tau_i(1);
tau_beta  = tau_i(2);

% Ricostruzione di Q_int tramite la relazione di Bryson
ximax = xmax_int .* tau_i;
Q_int = diag(1 ./ ximax.^2);

% Matrice di costo completa e sintesi LQR aumentata
Q_lqr = blkdiag(Q_heli, Q_act, Q_int);

K_aug = lqr(Aa, Ba, Q_lqr, R_lqr);

Krp = K_aug(:, 1:n_ext);
Kri = K_aug(:, n_ext+1:end);

Acl = Aa - Ba*K_aug;

% Verifiche
fprintf('\n=== TARATURA LQGI FINALE ===\n');
fprintf('tau_alpha = %.6f s\n', tau_alpha);
fprintf('tau_beta  = %.6f s\n', tau_beta);

disp('Q_int ottenuta:')
disp(Q_int)

disp('Poli del sistema aumentato:')
disp(eig(Acl))

assert(norm(diag(Q_int) - q_int_target, inf) < 1e-8, ...
    'I pesi Q_int non corrispondono ai valori desiderati.');

%% -------------------------------------------------
% 11) Costruzione controllore LQG 2DOF con integratore su alfa e beta
% --------------------------------------------------

ni = size(C_int,1);


Ac_ctrl = [
    zeros(ni,ni),              zeros(ni,n_ext);
    -B_ext*Kri,                A_ext-B_ext*Krp-Ke*C_ext
];

Bcr = [
    eye(ni);
    zeros(n_ext,ni)
];



Bcy = [
    -HalphaBeta;
    Ke
];


Cc_ctrl = [
    -Kri   -Krp;
    zeros(n_ext,ni)  eye(n_ext)
];

Dc_ctrl = [
    zeros(m,ni+c);
    zeros(n_ext,ni+c)
];

K_lqg_int = ss( ...
    Ac_ctrl, ...
    [Bcr Bcy], ...
    Cc_ctrl, ...
    Dc_ctrl);


K_lqg_int.InputName = {
    'r_alpha'
    'r_beta'   % Aggiunto per il secondo integratore
    'y_acc'
    'mx'
    'my'
};
K_lqg_int.OutputName = {
    'u1'
    'u2'
    'alpha_hat'
    'alphadot_hat'
    'beta_hat'
    'betadot_hat'
    'act1_hat'
    'act1dot_hat'
    'act1_pade_hat'
    'act2_hat'
    'act2dot_hat'
    'act2_pade_hat'
};

save('LQG_Controllers.mat', ...
    'K_lqg_int', 'Q_int', 'tau_i', '-append');

disp('==============================================')
disp(' K_lqg_int creato correttamente')
disp(' Inserire il blocco in Simulink')
size(K_lqg_int)
disp('==============================================')

%% ==========================================
% Connessione LQG 2DOF con riferimento
% ==========================================


%% Connessione nominale

Pnom = P_ext;

Pnom.InputName = {'u1','u2'};
Pnom.OutputName = {
    'y_acc'
    'mx'
    'my'
};

K_ctrl = K_lqg_int(1:2,:);

K_ctrl.InputName = {
    'r_alpha'
    'r_beta'
    'y_acc'
    'mx'
    'my'
};

K_ctrl.OutputName = {
    'u1'
    'u2'
};


CL_nom = connect(Pnom,K_ctrl,...
    {'r_alpha','r_beta'},...
    {'y_acc','mx','my'});

C_monitor = [
    C_ext;
    1 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0
];

D_monitor = zeros(5,2);

P_monitor = ss(A_ext,B_ext,C_monitor,D_monitor);

P_monitor.InputName = {'u1','u2'};

P_monitor.OutputName = { 
    'y_acc'
    'mx'
    'my'
    'delta_alpha'
    'delta_beta'
};

CL_monitor = connect(P_monitor,...
                     K_ctrl,...
                     {'r_alpha','r_beta'},...
                     {'delta_alpha','delta_beta'});
dcgain(CL_nom)
figure
step(CL_nom)

grid on
title('Risposta sensori nominale')


%% ==========================================
% Analisi stabilità
% ==========================================

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
title('Poli e zeri Closed Loop LQG')


% Risposta al gradino
figure
step(CL_nom)
grid on
title('risposta sensori nominale')



%% Analisi prestazioni: r_alpha = deg2rad(3) rad, r_beta = 0 rad

t = linspace(0,20,2001);

r = zeros(length(t),2);
r(:,1) = deg2rad(3);     % r_alpha = deg2rad(3) rad
r(:,2) = 0;              % r_beta  = 0 rad

[y_mon,t_mon] = lsim(CL_monitor,r,t);

delta_alpha_resp = y_mon(:,1);
delta_beta_resp  = y_mon(:,2);

figure
plot(t_mon,delta_alpha_resp,'LineWidth',1.5)
hold on
plot(t_mon,delta_beta_resp,'LineWidth',1.5)
grid on
legend('\alpha','\beta')
xlabel('Tempo [s]')
ylabel('Angolo [rad]')
title('Tracking LQG: r_\alpha = deg2rad(3) rad, r_\beta = 0')
info_alpha = stepinfo(delta_alpha_resp,t_mon);

disp(info_alpha)

disp('Tempo di assestamento alpha:')
disp(info_alpha.SettlingTime)

disp('Overshoot alpha:')
disp(info_alpha.Overshoot)

% Per beta il riferimento è zero:
beta_final = delta_beta_resp(end);
beta_max = max(abs(delta_beta_resp));
beta_rms = rms(delta_beta_resp);

fprintf('\nPrestazioni beta (r_beta = 0):\n')
fprintf('Valore finale beta = %.6g rad\n',beta_final)
fprintf('Massimo |beta|      = %.6g rad\n',beta_max)
fprintf('RMS beta            = %.6g rad\n',beta_rms)


size(Aa)
size(Ba)
size(K_aug)
eig(Aa-Ba*K_aug)
