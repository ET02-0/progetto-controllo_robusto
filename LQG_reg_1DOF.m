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

% Estrazione dinamica attuatori (secondo ordine, senza Padé per la sintesi)
[A_act, B_act, C_act, D_act] = ssdata(G_act_2nd);
Actuator = ss(A_act, B_act, C_act, D_act);

% MIMO attuatori indipendenti (2 ingressi, 2 uscite)
Actuators_MIMO = blkdiag(Actuator, Actuator);

%% 2) Connessione per il Modello Esteso (8 Stati)
% P_sintesi = Impianto Elicottero * Attuatori MIMO
P_sintesi = P_nom * Actuators_MIMO;
[A_ext, B_ext, C_ext, D_ext] = ssdata(P_sintesi);

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
Q_heli = diag([1000 40 1500 40]);
Q_act  = diag([1, 0.1, 1, 0.1]);     % Pesi moderati sugli attuatori
Q_lqr  = blkdiag(Q_heli, Q_act);
R_lqr = 0.1*eye(2)

K_lqr = lqr(A_ext, B_ext, Q_lqr, R_lqr);
disp('Guadagni LQR calcolati con successo.');

%% 5) Sintesi Filtro di Kalman (Stima dello stato)
% Il rumore di processo entra attraverso gli ingressi (es. disturbo sui rotori/vento)
Gk  = B_ext;
W = diag([
    0.001
    0.001
]);

V = sensor.R;

Qn = W;
Rn = V;

Vinv = diag(1./diag(V));

C_angle_meas = C_ext(:,[1 3]);

HalphaBeta = ...
    (C_angle_meas.'*Vinv*C_angle_meas) \ ...
    (C_angle_meas.'*Vinv);
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


K_lqg_reg = ss(Ac_ctrl,...
               Bc_ctrl,...
               Cc_ctrl,...
               Dc_ctrl);
save('LQG_Controllers.mat','K_lqg_reg');

K_lqg_reg.InputName={'y_acc','m_x','m_y'};
K_lqg_reg.OutputName = {'u1','u2'};


disp('Controllore LQG regolatore 1-DOF creato');
Acl_reg = A_ext - B_ext*K_lqr;
Aobs = A_ext - Ke*C_ext;

disp('Poli LQR:')
disp(eig(Acl_reg))

disp('Poli osservatore:')
disp(eig(Aobs))