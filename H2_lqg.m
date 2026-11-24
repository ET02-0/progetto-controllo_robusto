%% =========================================================================
%  H2_LQG_FORM - Formulazione H2 equivalente all'LQG a 2 DOF 
%  Modello esteso a 10 stati: elicottero + attuatori con Pade'.
%  Ingressi esogeni  w = [eta_w; eta_v]  (rumori bianchi unitari)
%  Uscite da pesare  z = [Q^(1/2) x ; R^(1/2) u]
% =========================================================================
clc; close all;

req = {'A_ext','B_ext','C_ext','Q_lqr','R_lqr','W','V'};
for k = 1:numel(req)
    assert(exist(req{k},'var') == 1, ...
        'Variabile mancante: %s (esegui prima dataset e script LQG a 2 DOF)', req{k});
end

%% 1. DATI: gli stessi dell'LQG a 2 DOF
A  = A_ext;   B = B_ext;   C = C_ext;
nx = size(A,1);  nu = size(B,2);  ny = size(C,1);

Gn = B_ext;                      % rumore di processo agli ingressi degli attuatori 
nw = size(Gn,2);

Q  = Q_lqr(1:nx,1:nx);           
R  = R_lqr;

Q  = 0.5*(Q+Q.');   R  = 0.5*(R+R.');
Wc = 0.5*(W+W.');   Vc = 0.5*(V+V.');

fprintf('autovalori di Q: min = %.3e, max = %.3e\n', min(eig(Q)), max(eig(Q)));
assert(min(eig(Q)) > -1e-9*max(1,max(eig(Q))), 'Q ha autovalori negativi: controllare Q_lqr nel workspace.');
[Uq,Dq] = eig(Q);
Qh = sqrt(max(diag(Dq),0)) .* Uq.';   % Qh'*Qh = Q anche se Q e' solo semidefinita
Rh = chol(R);                    % Rh'*Rh = R
Wh = chol(Wc,'lower');           % Wh*Wh' = W
Vh = chol(Vc,'lower');           % Vh*Vh' = V

%% 2. GENERALIZED PLANT H2 
B1  = [Gn*Wh, zeros(nx,ny)];
C1  = [Qh; zeros(nu,nx)];
D12 = [zeros(nx,nu); Rh];
D21 = [zeros(ny,nw), Vh];

P_H2 = ss(A, [B1 B], [C1; C], ...
          [zeros(nx+nu,nw+ny), D12; D21, zeros(ny,nu)]);

%% 3. SINTESI H2  (convenzione h2syn: u = K*y)
[K_H2_y, ~, gam_H2] = h2syn(P_H2, ny, nu);
K_H2_y = ss(K_H2_y);

fprintf('\n=============================================\n');
fprintf(' H2\n');
fprintf('=============================================\n');
fprintf('Ordine K_H2_y  = %d\n', order(K_H2_y));
fprintf('||Tzw||_2 H2   = %.10g\n', gam_H2);

%% 4. VERIFICA LQG
Klqr = lqr(A, B, Q, R);
Lkal = lqe(A, Gn, C, Wc, Vc);
K_LQG_y = ss(A - B*Klqr - Lkal*C, Lkal, -Klqr, zeros(nu,ny));

CL_H2  = lft(P_H2, K_H2_y);
CL_LQG = lft(P_H2, K_LQG_y);

gam_LQG = norm(CL_LQG, 2);
relGap  = abs(gam_H2 - gam_LQG) / max(gam_H2, eps);

fprintf('\nVERIFICA EQUIVALENZA (mappa w -> z)\n');
fprintf('||Tzw||_2 H2   = %.10g\n', gam_H2);
fprintf('||Tzw||_2 LQG  = %.10g\n', gam_LQG);
fprintf('Scarto relativo = %.3e\n', relGap);
if relGap < 1e-4
    disp('Equivalenza LQG-H2 verificata.');
else
    warning('Scarto non trascurabile: controllare convenzioni di segno e dati.');
end

figure('Name','Equivalenza LQG - H2');
sigma(CL_H2, CL_LQG, logspace(-3,3,600)); grid on;
legend('H_2 via h2syn','LQG (lqr+lqe)','Location','best');
title('Mappa w \rightarrow z: H_2 vs LQG');

%% 5. CONTROLLORE 2 DOF
%  u  = Kr*r + Ck*xk
%  xk' = Ak*xk + B*Kr*r + Bk*y
iA = find(any(abs(C) > 1e-12, 1));              
assert(numel(iA) == 2, 'C_ext dovrebbe avere due sole colonne non nulle (alpha, beta).');
C_ang = zeros(2,nx);
C_ang(1,iA(1)) = 1;
C_ang(2,iA(2)) = 1;

[Ak,Bk,Ck,Dk] = ssdata(K_H2_y);
coordOK = (norm(Dk) < 1e-9) && ...
          (norm(Ak - (A + B*Ck - Bk*C)) < 1e-6*max(1,norm(Ak)));
fprintf('\nStato di K_H2_y = stima dello stato dell''impianto: %d\n', coordOK);

if ~coordOK
    warning(['Realizzazione di h2syn non in coordinate di stato: uso quella ', ...
             'di K_LQG_y (stessa funzione di trasferimento, verificata sopra).']);
    [Ak,Bk,Ck,Dk] = ssdata(K_LQG_y);
end

Kr = -inv(C_ang*((A + B*Ck)\B));                
fprintf('Precompensatore Kr =\n'); disp(Kr);

nx_k = size(Ak, 1);  

% Amplio C aggiungendo l'identità 10x10 in basso (diventa 12x10)
Ck_ext = [Ck; eye(nx_k)]; 

% Amplio D aggiungendo righe di zeri in basso (diventa 12x5)
% Gli ingressi sono 5 (2 riferimenti + ny misure, con ny=3)
Dk_ext = [ [Kr, zeros(nu, ny)]; zeros(nx_k, 2 + ny) ];

% Assegno le nuove matrici Ck_ext e Dk_ext
K_H2_lqg = ss(Ak, [B*Kr, Bk], Ck_ext, Dk_ext);
K_H2_lqg.InputName  = {'r_alpha','r_beta','y_acc','m_x','m_y'};

% Aggiorno i nomi delle uscite per - "connect"
out_names = {'u1','u2'};
for i = 1:nx_k
    out_names{end+1} = sprintf('xk_%d', i); 
end
K_H2_lqg.OutputName = out_names;

% Anello chiuso nominale: r -> [alpha, beta]
Gall = ss(A, B, [C_ang; C], zeros(2+ny,nu));
Gall.InputName  = {'u1','u2'};
Gall.OutputName = {'alpha','beta','y_acc','m_x','m_y'};

CL_r = connect(Gall, K_H2_lqg, {'r_alpha','r_beta'}, {'alpha','beta'});
assert(isstable(CL_r), 'Anello chiuso nominale instabile.');
fprintf('Guadagni statici r -> angoli (devono essere ~ I):\n'); disp(dcgain(CL_r));

% Confronto diretto con il controllore LQG a 2 DOF, se presente
if exist('K_lqg_2dof','var')
    Kref   = K_lqg_2dof(1:2,:);
    CL_ref = connect(Gall, Kref, {'r_alpha','r_beta'}, {'alpha','beta'});
    fprintf('||CL_H2 - CL_LQG2dof||_inf = %.3e\n', norm(CL_r - CL_ref, inf));
end

t = 0:0.01:20;
figure('Name','H2-LQG: gradino di riferimento');
step(deg2rad(3)*CL_r, t); grid on;
title('H_2 (forma LQG): risposta a gradino di 3 deg (r \rightarrow \alpha,\beta)');

%% 6. VARIANTE CON MISURE ANGOLARI: stesso costo H2
% La ricostruzione WLS Hab*y_imu e' una statistica sufficiente: il costo H2
% ottimo non cambia se al controllore si danno gli angoli ricostruiti.
Cm   = C(:,iA);
Hab  = (Cm.'*(Vc\Cm)) \ (Cm.'/Vc);              
Vang = Hab*Vc*Hab.';  Vang = 0.5*(Vang+Vang.');
Vah  = chol(Vang,'lower');

P_ang = ss(A, [[Gn*Wh, zeros(nx,2)], B], [C1; C_ang], ...
           [zeros(nx+nu,nw+2), D12; [zeros(2,nw), Vah], zeros(2,nu)]);
[~, ~, gam_ang] = h2syn(P_ang, 2, nu);

fprintf('\nVARIANTE ANGOLARE\n');
fprintf('||Tzw||_2 (angoli) = %.10g  (scarto rel. vs IMU = %.2e)\n', ...
    gam_ang, abs(gam_ang-gam_H2)/gam_H2);

%% 7. SALVATAGGIO
save('H2_LQG_controller.mat', ...
     'K_H2_lqg','K_H2_y','K_LQG_y','Kr','gam_H2','gam_LQG','gam_ang', ...
     'relGap','P_H2','Hab');
disp('Salvato H2_LQG_controller.mat');