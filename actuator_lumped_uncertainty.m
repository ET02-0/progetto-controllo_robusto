%% ========================================================================
% INCERTEZZA MOLTIPLICATIVA CONCENTRATA DEGLI ATTUATORI
% Elicottero 2-DoF - adattato a dataset_elicottero
%
% Dall'incertezza parametrica dell'attuatore:
%
%       G_act_unc
%
% si costruisce una rappresentazione moltiplicativa:
%
%       G_act_unc_lumped = G_act_nom * (I + WI*Delta)
%
% con:
%       ||Delta||_inf <= 1
%
% Tale rappresentazione viene utilizzata per la mu-synthesis.
% ========================================================================

close all;
clc;

fprintf('============================================================\n');
fprintf(' INCERTEZZA MOLTIPLICATIVA CONCENTRATA ATTUATORI\n');
fprintf('============================================================\n');

%% ========================================================================
% 0. CARICAMENTO DEI MODELLI E DELLA NORMALIZZAZIONE
% ========================================================================

% Il dataset_elicottero deve essere stato eseguito prima.
requiredVariables = {
    'G_act_nom'
    'G_actuator_unc'
    'P_unc'
};

for k = 1:numel(requiredVariables)

    if ~exist(requiredVariables{k},'var')

        error( ...
            'Variabile mancante: %s. Eseguire prima dataset_elicottero.', ...
            requiredVariables{k});

    end

end

% Du, Dy, Du_inv e Dy_inv sono definiti in HINF_SETUP
% e salvati in HINF_workspace.mat.

if ~exist('HINF_workspace.mat','file')
    error('File HINF_workspace.mat non trovato. Eseguire prima HINF_SETUP.');
end

load('HINF_workspace.mat', ...
     'Du', ...
     'Dy', ...
     'Du_inv', ...
     'Dy_inv');

rng(10);
%% ========================================================================
% 1. ATTUATORI NOMINALI E INCERTI
% ========================================================================

G1_nom = G_act_nom;
G2_nom = G_act_nom;

G1_unc_param = G_actuator_unc;
G2_unc_param = G_actuator_unc;

%% ========================================================================
% 2. CAMPIONAMENTO DELL'INCERTEZZA PARAMETRICA
% ========================================================================

Nsample = 200;

G1_samples = usample( ...
    G1_unc_param, ...
    Nsample);

G2_samples = usample( ...
    G2_unc_param, ...
    Nsample);

fprintf('\nCampionati generati per ciascun attuatore: %d\n',Nsample);

%% ========================================================================
% 3. FIT DELL'INCERTEZZA MOLTIPLICATIVA
% ========================================================================

% Ordine del peso dinamico che ricopre l'errore relativo.

OrderWt = 2;

[~,info1] = ucover( ...
    G1_samples, ...
    G1_nom, ...
    OrderWt, ...
    'InputMult');

[~,info2] = ucover( ...
    G2_samples, ...
    G2_nom, ...
    OrderWt, ...
    'InputMult');

WI1 = minreal( ...
    tf(info1.W1), ...
    1e-7);

WI2 = minreal( ...
    tf(info2.W1), ...
    1e-7);

fprintf('\n============================================================\n');
fprintf(' PESI DI INCERTEZZA MOLTIPLICATIVA\n');
fprintf('============================================================\n');

disp('WI1 =');
WI1

disp('WI2 =');
WI2

%% ========================================================================
% 4. VERIFICA DEL COVER
% ========================================================================

omega = logspace(-1,3,500);

Rel1 = ...
    (G1_samples - G1_nom) / G1_nom;

Rel2 = ...
    (G2_samples - G2_nom) / G2_nom;

figure('Name','Attuatore 1 - Incertezza moltiplicativa');
sigma(Rel1,omega);
hold on;
sigma(WI1,omega);
grid on;
title( ...
    'Attuatore 1: errore relativo e peso W_{I,1}(j\omega)', ...
    'Interpreter','tex');

figure('Name','Attuatore 2 - Incertezza moltiplicativa');
sigma(Rel2,omega);
hold on;
sigma(WI2,omega);
grid on;
title( ...
    'Attuatore 2: errore relativo e peso W_{I,2}(j\omega)', ...
    'Interpreter','tex');

%% ========================================================================
% 5. COSTRUZIONE DEI BLOCCHI ULTIDYN
% ========================================================================

Delta_act1 = ultidyn( ...
    'Delta_act1', ...
    [1 1]);

Delta_act2 = ultidyn( ...
    'Delta_act2', ...
    [1 1]);

%% ========================================================================
% 6. MODELLO LUMPED DEGLI ATTUATORI
% ========================================================================

G1_unc_lumped = ...
    G1_nom * ...
    (1 + WI1*Delta_act1);

G2_unc_lumped = ...
    G2_nom * ...
    (1 + WI2*Delta_act2);

Gact_lumped = blkdiag( ...
    G1_unc_lumped, ...
    G2_unc_lumped);

Gact_lumped.InputName = {
    'u1'
    'u2'
};

Gact_lumped.OutputName = {
    'F1'
    'F2'
};

%% ========================================================================
% 7. RIDUZIONE DELLE INCERTEZZE MECCANICHE
% ========================================================================

% Per la mu-synthesis manteniamo solo le incertezze meccaniche
% che scegliamo di includere nel problema finale.
%
% N.B. Queste sono esattamente le ureal definite nel dataset_elicottero.

Pmech_reduced = P_unc;

parametersToNominal = {
    'Jy'
    'Jz'
    'm'
    'l'
};

for k = 1:numel(parametersToNominal)

    parName = parametersToNominal{k};

    if isfield( ...
            Pmech_reduced.Uncertainty, ...
            parName)

        block = ...
            Pmech_reduced.Uncertainty.(parName);

        Pmech_reduced = ...
            usubs( ...
                Pmech_reduced, ...
                parName, ...
                block.NominalValue);

    end

end

%% ========================================================================
% 8. ESTRAZIONE DELLE USCITE ANGOLARI
% ========================================================================

% P_unc è il modello con uscite:
%
%   delta_y_acc
%   delta_mx
%   delta_my
%
% Per la mu-synthesis ci servono soltanto:
%
%   delta_alpha
%   delta_beta
%
% quindi costruiamo direttamente il modello angolare.

C_angles = [
    1 0 0 0;
    0 0 1 0
];

Pq_mech_reduced = ss( ...
    Pmech_reduced.A, ...
    Pmech_reduced.B, ...
    C_angles, ...
    zeros(2,2));

Pq_mech_reduced.InputName = {
    'F1'
    'F2'
};

Pq_mech_reduced.OutputName = {
    'alpha'
    'beta'
};

%% ========================================================================
% 9. MODELLO FINALE PER MU-SYNTHESIS
% ========================================================================

G_uncertain_lumped = ...
    Pq_mech_reduced * ...
    Gact_lumped;

G_uncertain_lumped.InputName = {
    'u1'
    'u2'
};

G_uncertain_lumped.OutputName = {
    'alpha'
    'beta'
};

%% ========================================================================
% 10. NORMALIZZAZIONE
% ========================================================================

% Usiamo la stessa normalizzazione già definita nel dataset/setup:
%
%       G_scaled = Dy^-1 * G * Du

G_uncertain_lumped_scaled = ...
    Dy_inv * ...
    G_uncertain_lumped * ...
    Du;

G_uncertain_lumped_scaled.InputName = {
    'ubar1'
    'ubar2'
};

G_uncertain_lumped_scaled.OutputName = {
    'alpha_bar'
    'beta_bar'
};

%% ========================================================================
% 11. VERIFICA
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' MODELLO FINALE MU\n');
fprintf('============================================================\n');

disp('G_uncertain_lumped =');
G_uncertain_lumped

disp('G_uncertain_lumped_scaled =');
G_uncertain_lumped_scaled

%% ========================================================================
% 12. SALVATAGGIO
% ========================================================================

save( ...
    'ACTUATOR_LUMPED.mat', ...
    'WI1', ...
    'WI2', ...
    'Delta_act1', ...
    'Delta_act2', ...
    'Gact_lumped', ...
    'Pmech_reduced', ...
    'Pq_mech_reduced', ...
    'G_uncertain_lumped', ...
    'G_uncertain_lumped_scaled');

fprintf('\n============================================================\n');
fprintf(' ACTUATOR_LUMPED.mat SALVATO CORRETTAMENTE\n');
fprintf('============================================================\n');