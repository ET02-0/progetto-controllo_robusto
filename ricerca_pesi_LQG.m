%% Pesi iniziali da regola di Bryson

% Limiti di progetto degli stati:
% [alpha; alphadot; beta; betadot]
% Angoli in rad, velocita' angolari in rad/s.

xmax = [deg2rad(10);   % alpha: 10°
        0.75;           % alpha_dot: 0.75 rad/s
        deg2rad(7);   % beta: 7°
        0.5  ];          % beta_dot: 0.5 rad/s

Q_heli_base = diag(1 ./ xmax.^2);

rhoR =  1.7783;
R_fixed = rhoR * eye(2);

% Prima prova: non penalizziamo direttamente le coordinate
% interne degli attuatori, che non sono necessariamente
% grandezze fisiche normalizzate.
Q_act_fixed = zeros(6);

d0 = diag(Q_heli_base);

disp('Q_heli_base da Bryson:');
disp(Q_heli_base);

disp('Q_act_fixed:');
disp(Q_act_fixed);

disp('R_fixed = rhoR * eye(2), parametro di taratura:');
disp(R_fixed);

%% Griglia attorno alla base di Bryson

sa_values = logspace(-1,2,7);
sb_values = logspace(-1,1,5);

dAmp = [5e-3, 2e-3];

% Nessun integratore nella prima fase
Ki0 = zeros(size(B_ext,2),0);
%% Adattamento ai nomi del progetto
C = C_ext;
Vpsd = V;
%% Matrice dei disturbi aerodinamici
% Stati: [alpha, alphadot, beta, betadot, attuatori...]
% Disturbi: [tau_alpha, tau_beta]

Bd = zeros(size(A_ext,1),2);

Bd(2,1) = 1/J_alpha_nom;
Bd(4,2) = 1/Jbeta_nom;
% Verifiche dimensionali
assert(numel(d0) == 4, ...
    'Q_heli deve avere 4 pesi diagonali.');

assert(size(Q_act_fixed,1) == size(A_ext,1)-4 && ...
       size(Q_act_fixed,2) == size(A_ext,1)-4, ...
    'Q_act non è compatibile con il modello esteso.');

assert(size(B_ext,2) == 2 && isequal(size(R_fixed),[2 2]), ...
    'Sono attesi 2 ingressi e una R 2x2.');

assert(size(C,2) == size(A_ext,1) && ...
       size(Bd,1) == size(A_ext,1), ...
    'C o Bd non sono compatibili con A_ext.');

% Colonne:
% sa, sb, stable, peakA, peakB, rmsB, rmsU,
% peakU_normalizzato, u1Min, u1Max, u2Peak
results = NaN(numel(sa_values)*numel(sb_values),11);
k = 0;

for ia = 1:numel(sa_values)

    sa = sa_values(ia);

    for ib = 1:numel(sb_values)

        sb = sb_values(ib);
        k = k + 1;

        % Scala alpha e alpha_dot insieme,
        % beta e beta_dot insieme
        Qh = diag(d0(:).*[sa; sa; sb; sb]);

        % I pesi degli attuatori restano invariati
        Qx = blkdiag(Qh,Q_act_fixed);

        % Nuovo guadagno LQR
        Kx = lqr(A_ext,B_ext,Qx,R_fixed);

        % Valutazione analitica ad anello chiuso
        m = evalLQG(A_ext,B_ext,C,Bd,Ke, ...
                    Kx,Ki0,HalphaBeta,Vpsd,dAmp);

        results(k,:) = [sa, sb, double(m.stable), ...
                        m.peakA, m.peakB, ...
                        m.rmsB, m.rmsU, m.peakU, ...
                        m.u1Min, m.u1Max, m.u2Peak];

    end
end

Tgrid = array2table(results, ...
    'VariableNames',{'sa','sb','stable', ...
    'peakA_rad','peakB_rad','rmsB','rmsU', ...
    'peakU_norm','u1_min','u1_max','u2_absmax'});

Tgrid.stable = logical(Tgrid.stable);

Tgrid.peakA_deg = rad2deg(Tgrid.peakA_rad);
Tgrid.peakB_deg = rad2deg(Tgrid.peakB_rad);

Tgrid = sortrows(Tgrid,{'sb','sa'});

disp(Tgrid);

%% Riferimento per il vincolo sul rumore

Ki0 = zeros(size(B_ext,2),0);

Kx_base = lqr(A_ext, B_ext, ...
    blkdiag(Q_heli_base, Q_act_fixed), R_fixed);

m_base = evalLQG(A_ext, B_ext, C, Bd, Ke, ...
    Kx_base, Ki0, HalphaBeta, Vpsd, dAmp);


%% Criteri di ammissibilità

% Tolleranza massima del 20% rispetto alla base
rmsB_max = 1.2 * m_base.rmsB;
rmsU_max = 1.2 * m_base.rmsU;

ammissibili = Tgrid( ...
    Tgrid.stable & ...
    Tgrid.rmsB <= 1.2*m_base.rmsB & ...
    Tgrid.rmsU <= 1.2*m_base.rmsU, :);

ammissibili = sortrows(ammissibili, 'peakB_rad');

disp('Configurazioni ammissibili:');
disp(ammissibili);

%% Visualizza i 10 candidati migliori

nMostra = min(10, height(ammissibili));

disp('I migliori candidati ammissibili:');
disp(ammissibili(1:nMostra, ...
    {'sa','sb','peakA_deg','peakB_deg', ...
     'rmsB','rmsU','peakU_norm', ...
     'u1_min','u1_max','u2_absmax'}));

fprintf('Limite RMS beta  = %.6f\n', 1.2*m_base.rmsB);
fprintf('Limite RMS comando = %.6f\n', 1.2*m_base.rmsU);

%% Grafici per la scelta dei pesi

figure;

scatter(Tgrid.rmsB, Tgrid.peakB_deg, ...
        65, log10(Tgrid.sb), 'filled');

grid on;
xlabel('rmsB - norma H_2 rumore IMU \rightarrow beta');
ylabel('Picco beta [deg]');
title('Picco di beta rispetto al rumore in uscita');
cb = colorbar;
cb.Label.String = 'log10(sb)';

figure;

scatter(Tgrid.rmsU, Tgrid.peakB_deg, ...
        65, log10(Tgrid.sb), 'filled');

grid on;
xlabel('rmsU - norma H_2 rumore IMU \rightarrow comando');
ylabel('Picco beta [deg]');
title('Compromesso tra reiezione del disturbo e rumore sul comando');
cb = colorbar;
cb.Label.String = 'log10(sb)';
fprintf('\n=== VERIFICA PESI ===\n');
disp('File individuato da MATLAB:');
disp(which('ricerca_pesi_LQG'));

disp('Q_heli_base:');
disp(Q_heli_base);

disp('Q_act_fixed:');
disp(Q_act_fixed);

disp('R_fixed:');
disp(R_fixed);

disp('d0:');
disp(d0.');

%% VERIFICA INDIPENDENTE DEL CANDIDATO sa=1, sb=1

Qh_check = Q_heli_base;
Qx_check = blkdiag(Qh_check, Q_act_fixed);

Kx_check = lqr(A_ext, B_ext, Qx_check, R_fixed);

m_check = evalLQG(A_ext, B_ext, C_ext, Bd, Ke, ...
    Kx_check, zeros(size(B_ext,2),0), ...
    HalphaBeta, V, dAmp);

fprintf('\n===== VERIFICA INDIPENDENTE =====\n');

disp('Diagonale Qx del candidato base:');
disp(diag(Qx_check).');

disp('Matrice R utilizzata:');
disp(R_fixed);

fprintf('Picco alpha = %.6f gradi\n', rad2deg(m_check.peakA));
fprintf('Picco beta  = %.6f gradi\n', rad2deg(m_check.peakB));
fprintf('rmsB        = %.6f\n', m_check.rmsB);
fprintf('rmsU        = %.6f\n', m_check.rmsU);

disp('Riga corrispondente nella griglia:');

idx = abs(Tgrid.sa-1) < 1e-12 & ...
      abs(Tgrid.sb-1) < 1e-12;

disp(Tgrid(idx, ...
    {'sa','sb','peakA_deg','peakB_deg','rmsB','rmsU'}));

sa_star = ammissibili.sa(1);
sb_star = ammissibili.sb(1);

Qh_star = diag(d0(:) .* ...
    [sa_star; sa_star; sb_star; sb_star]);

Qx_star = blkdiag(Qh_star, Q_act_fixed);

Kx_star = lqr(A_ext, B_ext, Qx_star, R_fixed);

fprintf('Candidato scelto: sa = %.5g, sb = %.5g\n', ...
    sa_star, sb_star);

disp('Guadagno LQR Kx_star:');
disp(Kx_star);