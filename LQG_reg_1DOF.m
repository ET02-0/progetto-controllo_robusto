%% =========================================================
%  SINTESI LQG 1-DOF SENZA INTEGRATORE (8 STATI)
% =========================================================
% Nota: NON mettere clear all se vuoi mantenere il dataset attivo nel workspace!
close all; 

% Puliamo solo le variabili del controllore per evitare conflitti
clear Ac_ctrl Bc_ctrl Cc_ctrl Dc_ctrl K_lqg_2dof CL_nom CL_unc Kr;

disp('==============================================')
disp(' SINTESI LQG 1-DOF SENZA INTEGRATORE (8 STATI)')
disp('==============================================')
umax = 5; % esempio Nm
%% 1) Estrazione imanto nominale e attuatori dal Dataset
% P_nom è il modello nominale dell'elicottero:
% 4 stati [alpha, alpha_dot, beta, beta_dot],
% 2 ingressi di controllo e 3 uscite sensoriali
[A_heli, B_heli, C_heli, D_heli] = ssdata(P_nom);

n_ext = size(A_ext, 1); % Dovrebbe essere 8
m     = size(B_ext, 2); % 2
c     = size(C_ext, 1); % 2

fprintf('Stati modello esteso = %d (Attesi: 8)\n', n_ext);
if n_ext ~= 8
    error('Errore: il modello esteso non ha 8 stati.');
end

%% 3) Analisi di Controllabilità e Osservabilità
if rank(ctrb(A_ext, B_ext)) == n_ext
    disp('Il sistema esteso è CONTROLLABILE.');
else
    warning('Il sistema NON è completamente controllabile.');
end

if rank(obsv(A_ext, C_ext)) == n_ext
    disp('Il sistema esteso è OSSERVABILE.');
else
    warning('Il sistema NON è completamente osservabile.');
end

%% 4) Sintesi LQR (Regolazione)
% Stati: [alpha, alpha_dot, beta, beta_dot, act1, act1_dot, act2, act2_dot]
Q_heli = diag([800, 20, 10000, 500]);

Q_act = diag([5, 0.5, 5, 0.5]);

Q_lqr = blkdiag(Q_heli, Q_act);

R_lqr = diag([0.5 0.25]);

K_lqr = lqr(A_ext, B_ext, Q_lqr, R_lqr);
disp('Guadagni LQR calcolati con successo.');

%% 5) Sintesi Filtro di Kalman (Stima dello stato)
% Il rumore di processo entra attraverso gli ingressi (es. disturbo sui rotori/vento)
Gk = B_ext;

W = diag([0.003 0.0015]);

%{
% --- W sovrascrivibile dall'esterno (es. sweep_kalman_W.m) ---
% Se W esiste gia' nel base workspace (impostato dallo sweep), la si usa cosi' com'e'.
% Altrimenti si usa il valore di default per l'esecuzione singola dello script.
if ~exist('W','var')
    W = diag([0.01  0.01 ]);
end
%}


V = sensor.R;

Qn = W;
Rn = V;

fprintf('LQG2 -> W = [%.6g %.6g]\n', W(1,1), W(2,2));

Estimator = ss(A_ext, Gk, C_ext, zeros(c, m));
[~, Ke, ~] = kalman(Estimator, Qn, Rn);
disp('Filtro di Kalman calcolato con successo.');

%% =========================================================
% 6) Assemblaggio LQG REGOLATORE 1-DOF (Senza Integratore)
%
% Struttura:
%
%        y
%        |
%        v
%     Kalman
%        |
%       x_hat
%        |
%       -K
%        |
%        u
%
% Legge di controllo:
%
%       u = -K_lqr*x_hat
%
% Nessun riferimento.
% Usato per:
% - stabilizzazione equilibrio
% - reiezione disturbi
% - analisi robustezza
%
% =========================================================


Ac_ctrl = A_ext - B_ext*K_lqr - Ke*C_ext;

% ingresso = misura y
Bc_ctrl = Ke;

% uscita = comando u
Cc_ctrl = -K_lqr;

Dc_ctrl = zeros(m,c);


K_lqg_reg = ss( ...
    Ac_ctrl, ...
    Bc_ctrl, ...
    [Cc_ctrl;
     eye(n_ext)], ...
    [Dc_ctrl;
     zeros(n_ext,c)]);

K_lqg_reg.InputName = {'y_acc','m_x','m_y'};

K_lqg_reg.OutputName = {
    'u1'
    'u2'
    'alpha_hat'
    'alphadot_hat'
    'beta_hat'
    'betadot_hat'
    'act1_hat'
    'act1dot_hat'
    'act2_hat'
    'act2dot_hat'
};

save('LQG_Controllers.mat','K_lqg_reg');



disp('Controllore LQG regolatore 1-DOF creato');
Acl_reg = A_ext - B_ext*K_lqr;
Aobs = A_ext - Ke*C_ext;

disp('Poli LQR:')
disp(eig(Acl_reg))

disp('Poli osservatore:')
disp(eig(Aobs))

%% =========================================================
% GRAFICO POLI LQG 1-DOF
% =========================================================

figure('Name','Poli LQG 1-DOF','Color','w')

plot(real(eig(Acl_reg)),imag(eig(Acl_reg)),'x', ...
    'MarkerSize',10, ...
    'LineWidth',2)

hold on

plot(real(eig(Aobs)),imag(eig(Aobs)),'o', ...
    'MarkerSize',8, ...
    'LineWidth',1.5)

xline(0,'k--','LineWidth',0.8)

grid on

xlabel('Parte reale')
ylabel('Parte immaginaria')

title('Poli LQG 1-DOF')

legend( ...
    'Poli LQR', ...
    'Poli osservatore', ...
    'Asse immaginario', ...
    'Location','best')