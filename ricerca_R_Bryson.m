%% RICERCA R - BASE DI BRYSON

C = C_ext;
Vpsd = V;

%% Matrice dei disturbi aerodinamici

Bd = zeros(size(A_ext,1),2);

Bd(2,1) = 1/J_alpha_nom;
Bd(4,2) = 1/Jbeta_nom;

%% Base di Bryson per gli stati dell'elicottero

xmax = [deg2rad(3);
        deg2rad(10);
        deg2rad(1);
        deg2rad(10)];

Q_heli_base = diag(1 ./ xmax.^2);

% Per ora non penalizziamo direttamente le coordinate
% interne degli attuatori
deltaF1_max = 0.75;
deltaF2_max = 0.75;

% Matrice di uscita degli attuatori
[A_act, ~, C_act, D_act] = ssdata(Actuators_MIMO);
W_deltaF = diag(1 ./ [deltaF1_max; deltaF2_max].^2);
Q_act_fixed = C_act' * W_deltaF * C_act;

Q_base = blkdiag(Q_heli_base,Q_act_fixed);

%% Ricerca del fattore R

% Forma relativa degli ingressi fissata
R_shape = eye(2);

% Fattore comune da provare
rhoR_values = logspace(-1,1,9);

dAmp = [5e-3, 2e-3];
Ki0 = zeros(size(B_ext,2),0);

% Colonne:
% rhoR, stable, peakA_deg, peakB_deg,
% rmsA, rmsB, rmsU1, rmsU2,
% u1_min, u1_max, u2_absmax

Rresults = NaN(numel(rhoR_values),11);

for k = 1:numel(rhoR_values)

    rhoR = rhoR_values(k);
    R = rhoR * R_shape;

    Kx = lqr(A_ext,B_ext,Q_base,R);

    m = evalLQG(A_ext,B_ext,C,Bd,Ke, ...
                Kx,Ki0,HalphaBeta,Vpsd,dAmp);

    Rresults(k,:) = ...
        [rhoR, double(m.stable), ...
         rad2deg(m.peakA), rad2deg(m.peakB), ...
         m.rmsA, m.rmsB, m.rmsU1, m.rmsU2, ...
         m.u1Min, m.u1Max, m.u2Peak];

end

%% Tabella dei risultati

TR = array2table(Rresults, ...
    'VariableNames',{'rhoR','stable', ...
    'peakA_deg','peakB_deg','rmsA','rmsB', ...
    'rmsU1','rmsU2','u1_min','u1_max','u2_absmax'});

TR.stable = logical(TR.stable);

TR.rmsCommand = hypot(TR.rmsU1,TR.rmsU2);

disp('Ricerca di R con Q di Bryson fissata:');
disp(TR);

%% Grafico: picco beta e rumore sui comandi

figure;

scatter(TR.rmsCommand,TR.peakB_deg, ...
        75,log10(TR.rhoR),'filled');

grid on;
xlabel('Norma H_2 complessiva rumore IMU -> comandi');
ylabel('Picco beta [deg]');
title('Picco beta e rumore trasferito ai comandi');

cb = colorbar;
cb.Label.String = 'log10(rhoR)';

%% Grafico: picco beta e rumore su beta

figure;

scatter(TR.rmsB,TR.peakB_deg, ...
        75,log10(TR.rhoR),'filled');

grid on;
xlabel('Norma H_2 rumore IMU -> beta');
ylabel('Picco beta [deg]');
title('Picco beta e rumore in uscita');

cb = colorbar;
cb.Label.String = 'log10(rhoR)';