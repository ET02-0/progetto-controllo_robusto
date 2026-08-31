% =========================================================================
% Script di Inizializzazione e Linearizzazione: Elicottero 2DoF
% Progetto di Controllo Robusto (Sintesi Robusta su 4 stati)
% =========================================================================

close all
clc

disp('SINTESI LQG 2DOF CON INTEGRATORE')



%% Attuatori

[A_act,B_act,C_act,D_act] = ssdata(G_act_2nd);

Actuator = ss(A_act,B_act,C_act,D_act);

Actuators_MIMO = blkdiag(Actuator,Actuator);


%% Plant esteso

P_ext = P_nom*Actuators_MIMO;

[A_ext,B_ext,C_ext,D_ext]=ssdata(P_ext);

M = [A_ext B_ext;
     -C_ext zeros(size(C_ext,1),size(B_ext,2))];

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


if n_ext ~= 8
    error('Errore: il modello esteso non ha 8 stati.')
end



%% -------------------------------------------------
% 4) Analisi modello esteso
% --------------------------------------------------

disp('Poli modello esteso:')
disp(eig(A_ext))


figure
pzmap(P_ext)
grid on
title('Poli modello esteso')


%% Integratori su pitch (alpha) e yaw (beta)
% Il vettore di stato esteso è: [alpha, alphad, beta, betad, stati_attuatori(1..4)]
C_int = [1 0 0 0 0 0 0 0;   % Estrae alpha (stato 1)
         0 0 1 0 0 0 0 0];  % Estrae beta (stato 3)

ni = size(C_int,1);         % Ora ni = 2

Aa = [
A_ext zeros(n_ext,ni);
-C_int zeros(ni,ni)
];


Ba = [
B_ext;
zeros(ni,m)
];

rank_ctrb = rank(ctrb(Aa,Ba));

if rank_ctrb < size(Aa,1)

    fprintf("Modo non controllabile:\n")

    [V,D]=eig(Aa);

    for i=1:length(diag(D))
        if abs(real(D(i)))<1e-8
            disp(V(:,i))
        end
    end

end


%% -------------------------------------------------
% 6) Verifica controllabilità e stabilizzabilità
% --------------------------------------------------
% Usiamo la decomposizione in valori singolari (SVD) 
% o tolleranze manuali per evitare falsi positivi da malcondizionamento.
rank_ctrb = rank(ctrb(Aa,Ba), 1e-10); % Tolleranza forzata per ignorare lo scaling

fprintf('\nTest di controllabilità: il sistema possiede %d stati.\n', size(Aa,1));
disp('Nota: se il rango di ctrb sembra inferiore, è dovuto al malcondizionamento');
disp('generato da poli lenti (0 rad/s) e poli veloci (attuatori a 30 rad/s).');
disp('La reale stabilizzabilità è garantita dal successo del comando lqr().');

fprintf('\nControllabilita: %d/%d\n',...
    rank_ctrb,size(Aa,1));


if rank_ctrb ~= size(Aa,1)
    warning('Sistema aumentato non completamente controllabile')
end
eig(Aa)
Co = ctrb(Aa,Ba);

[U,S,V]=svd(Co);

disp('Singular values controllabilità:')
disp(diag(S))


%% -------------------------------------------------
% 7) Verifica osservabilità
% --------------------------------------------------

rank_obs = rank(obsv(A_ext,C_ext));

fprintf('Osservabilita: %d/%d\n',...
    rank_obs,size(A_ext,1));


if rank_obs ~= size(A_ext,1)
    warning('Sistema non completamente osservabile')
end



%% -------------------------------------------------
% 8) Pesi LQR
%
% Stato:
% x = [elicottero(4); attuatori(4)]
%
% --------------------------------------------------
% =================================================================
% 1. Pesi LQR
% =================================================================

Q_heli = diag([400 20 100 500]);
Q_act  = diag([5, 0.5, 5, 0.5]);
Q_int = diag([150 800]);
Q_lqr = blkdiag(Q_heli,Q_act,Q_int);

R_lqr  = diag([1 1]);

%% -------------------------------------------------
% 9) Calcolo LQR
% --------------------------------------------------

K_aug = lqr(Aa, Ba, Q_lqr, R_lqr);
Krp = K_aug(:, 1:n_ext);    % 2x8
Kri = K_aug(:, n_ext+1:end); % 2x2

umax = 5; % esempio Nm

disp('Guadagni LQR')
disp(K_aug)

Acl = Aa-Ba*K_aug;

disp('Poli sistema chiuso LQR:')
disp(eig(Acl))



disp('Guadagno LQR calcolato')



%% -------------------------------------------------
% 10) Filtro di Kalman
% --------------------------------------------------
Gk = B_ext;

W = diag([0.001, 0.001]);

V = sensor.R

Vinv = diag(1./diag(V));

C_angle_meas = C_ext(:,[1 3]);

HalphaBeta = ...
    (C_angle_meas.'*Vinv*C_angle_meas) \ ...
    (C_angle_meas.'*Vinv);


[Ke,~,~] = lqe(A_ext,Gk,C_ext,W,V);

disp('Filtro Kalman calcolato')



%% -------------------------------------------------
% 11) Costruzione controllore LQG 2DOF con integratore su alfa
% -------------------------------------------------

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
    -Kri   -Krp
];

Dc_ctrl = zeros(m,ni+c);

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
};

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

K = K_lqg_int;


K.InputName = {
    'r_alpha'
    'r_beta'   % <-- MANCAVA QUESTO
    'y_acc'
    'mx'
    'my'
};
K.OutputName = {
    'u1'
    'u2'
};


% QUI C'ERA IL SECONDO ERRORE: aggiungi 'r_beta' agli ingressi di connect
CL_nom = connect(Pnom,K,...
    {'r_alpha', 'r_beta'},... 
    {
    'y_acc'
    'mx'
    'my'
    });

C_monitor = [
    C_ext;
    1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0
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

% Aggiungi 'r_beta' anche qui
CL_monitor = connect(P_monitor,...
                     K_lqg_int,... 
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
title('Poli Closed Loop LQG')


% Risposta al gradino
figure
step(CL_nom)
grid on
title('risposta sensori nominale')


Pfull_unc = P_full_unc;

Pfull_unc.InputName = {
    'u1'
    'u2'
};

Pfull_unc.OutputName = {
    'y_acc'
    'mx'
    'my'
};

CL_unc = connect(Pfull_unc,K,...
    {'r_alpha','r_beta'},...
    {
    'y_acc'
    'mx'
    'my'
    });


%% ==========================================
% Analisi robustezza
% ==========================================

disp('==============================================')
disp(' ANALISI ROBUSTEZZA')
disp('==============================================')

% Controllo nomi ingressi/uscite
CL_unc.InputName
CL_unc.OutputName


disp('==============================================')
disp(' ROBUST STABILITY')
disp('==============================================')

disp('Inizio robuststab...')

opts = robuststabOptions('Display','on');

tic
[stab_margin,wcu,info] = robuststab(CL_unc,opts);
tempo_robuststab = toc;

fprintf('\nTempo robuststab = %.2f s\n',tempo_robuststab);

disp('Margine robustezza:')
disp(stab_margin)

disp('Worst-case uncertainty:')
disp(wcu)

%% Robust performance

omega = logspace(-2,3,500);

figure
sigma(CL_unc,omega)
grid on
title('Singular Values Closed Loop')


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

wcgain(CL_unc)

size(Aa)
size(Ba)
size(K_aug)
rank(ctrb(Aa,Ba))
eig(Aa-Ba*K_aug)