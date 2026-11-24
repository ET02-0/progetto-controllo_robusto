%% =========================================================================
%  H2_SOFT_INTEGRATOR - H2 con disturbance shaping a bassa frequenza
% =========================================================================
clc; close all;

req = {'A_ext','B_ext','C_ext','Q_lqr','R_lqr','W','V','J_alpha_nom','Jbeta_nom'};
for k = 1:numel(req)
    assert(exist(req{k},'var') == 1, ...
        'Variabile mancante: %s (esegui prima dataset e script LQG a 2 DOF)', req{k});
end

%% 1. DATI (come H2_LQG_form)
A  = A_ext;   B = B_ext;   C = C_ext;
nx = size(A,1);  nu = size(B,2);  ny = size(C,1);
Gn = B_ext;      nw = size(Gn,2);   nd = 2;

Q  = Q_lqr(1:nx,1:nx);   R = R_lqr;
Q  = 0.5*(Q+Q.');   R  = 0.5*(R+R.');
Wc = 0.5*(W+W.');   Vc = 0.5*(V+V.');
fprintf('autovalori di Q: min = %.3e, max = %.3e\n', min(eig(Q)), max(eig(Q)));
assert(min(eig(Q)) > -1e-9*max(1,max(eig(Q))), 'Q ha autovalori negativi: controllare Q_lqr nel workspace.');
[Uq,Dq] = eig(Q);
Qh = sqrt(max(diag(Dq),0)) .* Uq.';   % Qh'*Qh = Q anche se Q e' solo semidefinita
Rh = chol(R);
Wh = chol(Wc,'lower');   Vh = chol(Vc,'lower');

%% 2. INGRESSO DELLE COPPIE AERODINAMICHE SUGLI STATI
% alpha e beta sono le colonne non nulle di C_ext; ordine heli = [a, ad, b, bd]
iA = find(any(abs(C) > 1e-12, 1));
assert(numel(iA) == 2, 'C_ext dovrebbe avere due sole colonne non nulle (alpha, beta).');
C_ang = zeros(2,nx);
C_ang(1,iA(1)) = 1;
C_ang(2,iA(2)) = 1;

Bd = zeros(nx,nd);
Bd(iA(1)+1,1) = 1/J_alpha_nom;       % d_alpha -> alpha_ddot
Bd(iA(2)+1,2) = 1/Jbeta_nom;         % d_beta  -> beta_ddot

%% 3. PARAMETRI DEL SOFT INTEGRATOR
% Con xd_dot = -eps*xd + wd*eta, la deviazione standard stazionaria di xd e'
% wd/sqrt(2*eps). Si fissa quindi wd dalla std. desiderata della coppia.
soft.eps  = 1e-6;                    % [rad/s]
soft.rms  = [5e-3; 2e-3];            % [N*m] come le ampiezze aero.* del progetto
soft.gain = 250;                       % fattore di taratura (>1: piu' azione integrale)
soft.wd   = soft.gain * soft.rms * sqrt(2*soft.eps);

fprintf('\n=============================================\n');
fprintf(' H2 SOFT INTEGRATOR\n');
fprintf('=============================================\n');
fprintf('eps = %.2e rad/s, gain = %.2f\n', soft.eps, soft.gain);
fprintf('wd_alpha = %.3e, wd_beta = %.3e\n', soft.wd(1), soft.wd(2));

%% 4. GENERALIZED PLANT H2 AUMENTATO
% w = [eta_w; eta_v; eta_d]   z = [Q^(1/2) x ; R^(1/2) u]   y = C x + Vh*eta_v
Asf = [A, Bd; zeros(nd,nx), -soft.eps*eye(nd)];
B1  = [Gn*Wh,        zeros(nx,ny), zeros(nx,nd);
       zeros(nd,nw), zeros(nd,ny), diag(soft.wd)];
B2  = [B; zeros(nd,nu)];
C1  = [Qh, zeros(nx,nd); zeros(nu,nx+nd)];
C2  = [C, zeros(ny,nd)];
D12 = [zeros(nx,nu); Rh];
D21 = [zeros(ny,nw), Vh, zeros(ny,nd)];

P_soft = ss(Asf, [B1 B2], [C1; C2], ...
            [zeros(nx+nu,nw+ny+nd), D12; D21, zeros(ny,nu)]);

%% 5. SINTESI H2
[K_soft_y, ~, gam_soft] = h2syn(P_soft, ny, nu);
K_soft_y = ss(K_soft_y);

fprintf('Ordine K_soft_y = %d\n', order(K_soft_y));
fprintf('||Tzw||_2       = %.10g\n', gam_soft);
fprintf('Anello chiuso (generalized plant) stabile = %d\n', ...
    isstable(lft(P_soft, K_soft_y)));

%% 6. CONTROLLORE 2 DOF, STESSA STRUTTURA DI K_lqg_2dof
%  u  = Kr*r + Ck*xk
%  xk' = Ak*xk + B2*Kr*r + Bk*y       (r entra anche nello stimatore)
%  Lo stato xk e' la stima di [x; xd]: lo si verifica numericamente.
[Ak,Bk,Ck,Dk] = ssdata(K_soft_y);
coordOK = (norm(Dk) < 1e-9) && ...
          (norm(Ak - (Asf + B2*Ck - Bk*C2)) < 1e-6*max(1,norm(Ak)));
fprintf('\nStato di K_soft_y = stima di [x; xd]: %d\n', coordOK);

if ~coordOK
    warning(['Realizzazione di h2syn non in coordinate di stato: ', ...
             'ricostruisco lo stesso controllore con lqr/lqe sul modello aumentato.']);
    Fs = lqr(Asf, B2, blkdiag(Q, zeros(nd)), R);
    Ls = lqe(Asf, [Gn*Wh, zeros(nx,nd); zeros(nd,nw), diag(soft.wd)], ...
             C2, eye(nw+nd), Vc);
    Ak = Asf - B2*Fs - Ls*C2;   Bk = Ls;   Ck = -Fs;   Dk = zeros(nu,ny);
end

% Precompensatore: a regime e senza disturbo la stima di xd e' nulla
Kr = -inv(C_ang*((A + B*Ck(:,1:nx))\B));
fprintf('Precompensatore Kr =\n'); disp(Kr);

nx_k = size(Ak, 1); % È 12 nel tuo caso

% Amplio C per estrarre i primi 10 stati.
% eye(10, nx_k) crea una matrice 10x12 (estrae i primi 10 stati e ignora gli ultimi 2)
% Dimensioni Ck_ext: 12x12 (2 righe di Ck + 10 righe per le stime)
Ck_ext = [Ck; eye(10, nx_k)]; 

% Amplio D aggiungendo 10 righe di zeri in basso
% Dimensioni Dk_ext: 12x5
Dk_ext = [ [Kr, zeros(nu, ny)]; zeros(10, 2 + ny) ];

% 3. Creo il controllore
K_H2_soft = ss(Ak, [B2*Kr, Bk], Ck_ext, Dk_ext);
K_H2_soft.InputName  = {'r_alpha','r_beta','y_acc','m_x','m_y'};

% 4. Aggiorno i nomi delle uscite per avere u1, u2 e xk_1 ... xk_10
out_names = {'u1','u2'};
for i = 1:10
    out_names{end+1} = sprintf('xk_%d', i); 
end
K_H2_soft.OutputName = out_names;

Gall = ss(A, B, [C_ang; C], zeros(2+ny,nu));
Gall.InputName  = {'u1','u2'};
Gall.OutputName = {'alpha','beta','y_acc','m_x','m_y'};

CL_r = connect(Gall, K_H2_soft, {'r_alpha','r_beta'}, {'alpha','beta'});
assert(isstable(CL_r), 'Anello chiuso nominale instabile.');
fprintf('Guadagni statici r -> angoli (devono essere ~ I):\n'); disp(dcgain(CL_r));

%% 7. VERIFICHE NEL TEMPO
t = 0:0.01:100;

% (a) gradino di riferimento 3 deg
figure('Name','Soft integrator: riferimento');
step(deg2rad(3)*CL_r, t); grid on;
title('H_2 soft: risposta a gradino di 3 deg (r \rightarrow \alpha,\beta)');

% (b) gradino di coppia aerodinamica (ampiezze del progetto), r = 0
if exist('aero','var')
    Dd = diag([aero.alpha.amplitude, aero.beta.amplitude]);
else
    Dd = diag(soft.rms);
end
Sel   = [zeros(ny,2), eye(ny)];
Gd    = ss(A, [B Bd], [C_ang; C], zeros(2+ny,nu+nd));
Gcl_d = feedback(Gd, K_soft_y*Sel, 1:nu, 1:(2+ny), +1);
figure('Name','Soft integrator: disturbi');
step(Gcl_d(1:2, nu+1:nu+nd)*Dd, t); grid on;
title('H_2 soft: risposta ai disturbi di coppia d_\alpha, d_\beta');

% (c) confronto con LQG a 2 DOF sullo stesso disturbo, se presente
if exist('K_lqg_2dof','var')
    Kref   = K_lqg_2dof(1:2,:);
    Gd2    = ss(A, [B Bd], [C_ang; C], zeros(2+ny,nu+nd));
    Gd2.InputName  = {'u1','u2','d_alpha','d_beta'};
    Gd2.OutputName = {'alpha','beta','y_acc','m_x','m_y'};
    CLd_ref = connect(Gd2, Kref, {'r_alpha','r_beta','d_alpha','d_beta'}, {'alpha','beta'});
    figure('Name','Disturbo: LQG2dof vs H2 soft');
    step(CLd_ref(:,3:4)*Dd, Gcl_d(1:2,nu+1:nu+nd)*Dd, t); grid on;
    legend('LQG 2DOF','H_2 soft','Location','best');
    title('Risposta ai disturbi di coppia: LQG 2DOF vs H_2 soft');
end

% (d) forma spettrale del disturbo assunto
s = tf('s');
Wdist = blkdiag(soft.wd(1)/(s+soft.eps), soft.wd(2)/(s+soft.eps));
figure('Name','Soft integrator: disturbance shaping');
sigma(Wdist, logspace(-4,2,500)); grid on;
title('Disturbance shaping assunto (\eta_d \rightarrow d)');

%% 8. SALVATAGGIO
save('H2_soft_controller.mat', ...
     'K_H2_soft','K_soft_y','Kr','gam_soft','soft','P_soft','Bd');
disp('Salvato H2_soft_controller.mat');