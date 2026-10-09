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
     'Du_inv', ...
     'Dy_inv');

rng(10);

%% ========================================================================
% 1. ATTUATORI NOMINALI
% ========================================================================

G1_nom = G_act_nom;
G2_nom = G_act_nom;


%% ========================================================================
% 2. INVILUPPO DETERMINISTICO DELL'INCERTEZZA COMPLESSIVA
% ========================================================================

% Questa procedura comprende:
%   - omega_n variabile di +/-10%
%   - tau_d variabile di +/-20%
%   - intero contributo W_Pade*Delta_Pade
%
% Non utilizziamo usample per costruire il peso.

omega = logspace(-1,3,500);
jw = 1i*omega(:);

% Griglie dei parametri per la costruzione dell'inviluppo
wn_vec = linspace( ...
    0.90*omega_n, ...
    1.10*omega_n, ...
    17);

tau_vec = linspace( ...
    0.80*tau_d_nom, ...
    1.20*tau_d_nom, ...
    41);

% Risposte nominali e modulo del peso Padé
Gnom_resp = squeeze( ...
    freqresp(G_act_nom,omega));
Gnom_resp = Gnom_resp(:);

Wp_abs = squeeze( ...
    abs(freqresp(W_Pade,omega)));
Wp_abs = Wp_abs(:);

% Inviluppo dell'errore relativo comprensivo di Delta_Pade
Ngrid = numel(wn_vec)*numel(tau_vec);
E_grid = zeros(numel(omega),Ngrid);

idx = 0;

for iwn = 1:numel(wn_vec)

    wn = wn_vec(iwn);

    G2nd_resp = ...
        wn^2 ./ ...
        (jw.^2 + 2*zeta*wn*jw + wn^2);

    for itau = 1:numel(tau_vec)

        tau = tau_vec(itau);

        % Padé del primo ordine per il parametro tau
        Gpade_resp = ...
            (1 - tau*jw/2) ./ ...
            (1 + tau*jw/2);

        % Rapporto tra attuatore parametrico e nominale
        R = ...
            (G2nd_resp.*Gpade_resp) ./ Gnom_resp;

        idx = idx + 1;

        % Bound:
        % |R*(1+W_Pade*Delta_Pade)-1|
        % <= |R-1| + |R|*|W_Pade|
        E_grid(:,idx) = ...
            abs(R-1) + abs(R).*Wp_abs;

    end

end

% Massimo sulle combinazioni dei parametri
E_env = max(E_grid,[],2);

fprintf('\n============================================================\n');
fprintf(' INVILUPPO DETERMINISTICO ATTUATORI\n');
fprintf('============================================================\n');

fprintf('Combinazioni parametriche analizzate: %d\n',Ngrid);

fprintf('Massimo inviluppo = %.6f\n',max(E_env));


%% ========================================================================
% 3. FIT DEL PESO MOLTIPLICATIVO COMUNE
% ========================================================================

% Margine per il fitting e per la discretizzazione delle griglie.
% Il risultato sarà comunque sottoposto a validazione indipendente.

Margine_WI = 1.10;

E_fit = max( ...
    Margine_WI*E_env, ...
    1e-8);

E_fit_frd = frd( ...
    reshape(E_fit,1,1,[]), ...
    omega);

OrderWt = 4;

% Vincolo esplicito: |WI| >= inviluppo maggiorato
Cfit = struct();
Cfit.LowerBound = E_fit_frd;
Cfit.UpperBound = Inf;

WI_fit = fitmagfrd( ...
    E_fit_frd, ...
    OrderWt, ...
    [], ...
    [], ...
    Cfit);

WI = minreal( ...
    tf(WI_fit), ...
    1e-7);

% Stessa famiglia fisica per i due attuatori.
% Le incertezze Delta finali resteranno indipendenti.
WI1 = WI;
WI2 = WI;

fprintf('\n============================================================\n');
fprintf(' PESO DI INCERTEZZA LUMPED\n');
fprintf('============================================================\n');

disp('WI1 = WI');
WI1

disp('WI2 = WI');
WI2


%% ========================================================================
% 4. VERIFICA DEL COVER SULLA GRIGLIA DI FIT
% ========================================================================

mag_WI = squeeze( ...
    abs(freqresp(WI,omega)));
mag_WI = mag_WI(:);

ratio_fit = E_env./mag_WI;

fprintf('\n============================================================\n');
fprintf(' VERIFICA COVER SULLA GRIGLIA DI FIT\n');
fprintf('============================================================\n');

fprintf('Massimo rapporto inviluppo/WI = %.6f\n', ...
    max(ratio_fit));

fprintf('Frequenze oltre il cover = %d / %d\n', ...
    sum(ratio_fit > 1),numel(omega));

figure('Name','Cover deterministico attuatori');

semilogx( ...
    omega, ...
    E_env, ...
    'LineWidth',1.5);

hold on;

semilogx( ...
    omega, ...
    mag_WI, ...
    'LineWidth',1.5);

grid on;

xlabel('\omega [rad/s]');
ylabel('Modulo');

title('Inviluppo completo dell''incertezza e peso lumped');

legend( ...
    'Inviluppo errore relativo', ...
    '|W_I|', ...
    'Location','best');


%% ========================================================================
% 4.1 VALIDAZIONE SU GRIGLIA PIU' DENSA DEI PARAMETRI
% ========================================================================

omega_val = logspace(-1,3,1000);
jw_val = 1i*omega_val(:);

wn_val = linspace( ...
    0.90*omega_n, ...
    1.10*omega_n, ...
    31);

tau_val = linspace( ...
    0.80*tau_d_nom, ...
    1.20*tau_d_nom, ...
    61);

Gnom_val = squeeze( ...
    freqresp(G_act_nom,omega_val));
Gnom_val = Gnom_val(:);

Wp_abs_val = squeeze( ...
    abs(freqresp(W_Pade,omega_val)));
Wp_abs_val = Wp_abs_val(:);

mag_WI_val = squeeze( ...
    abs(freqresp(WI,omega_val)));
mag_WI_val = mag_WI_val(:);

E_env_val = zeros(numel(omega_val),1);

for iwn = 1:numel(wn_val)

    wn = wn_val(iwn);

    G2nd_resp = ...
        wn^2 ./ ...
        (jw_val.^2 + 2*zeta*wn*jw_val + wn^2);

    for itau = 1:numel(tau_val)

        tau = tau_val(itau);

        Gpade_resp = ...
            (1 - tau*jw_val/2) ./ ...
            (1 + tau*jw_val/2);

        R = ...
            (G2nd_resp.*Gpade_resp) ./ Gnom_val;

        E_candidate = ...
            abs(R-1) + abs(R).*Wp_abs_val;

        E_env_val = max( ...
            E_env_val, ...
            E_candidate);

    end

end

ratio_val = E_env_val./mag_WI_val;

fprintf('\n============================================================\n');
fprintf(' VALIDAZIONE SU GRIGLIA DENSA\n');
fprintf('============================================================\n');

fprintf('Coppie parametriche di validazione: %d\n', ...
    numel(wn_val)*numel(tau_val));

fprintf('Massimo rapporto inviluppo/WI = %.6f\n', ...
    max(ratio_val));

fprintf('Frequenze oltre il cover = %d / %d\n', ...
    sum(ratio_val > 1),numel(omega_val));

figure('Name','Validazione cover lumped');

semilogx( ...
    omega_val, ...
    E_env_val, ...
    'LineWidth',1.5);

hold on;

semilogx( ...
    omega_val, ...
    mag_WI_val, ...
    'LineWidth',1.5);

grid on;

xlabel('\omega [rad/s]');
ylabel('Modulo');

title('Validazione indipendente del peso lumped');

legend( ...
    'Inviluppo di validazione', ...
    '|W_I|', ...
    'Location','best');

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
%
% Per la mu-synthesis manteniamo soltanto:
%
%   J_alpha
%   l
%
% Le altre incertezze meccaniche vengono fissate al valore nominale:
%
%   Jy
%   Jz
%   m
%   eps_p
%   eps_y
%
% ================================================================

Pmech_reduced = P_unc;

parametersToNominal = {
    'Jy'
    'Jz'
    'm'
    'eps_p'
    'eps_y'
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

fprintf('\n============================================================\n');
fprintf('INCERTEZZE DOPO LA RIDUZIONE MECCANICA\n');
fprintf('============================================================\n');

disp(Pmech_reduced.Uncertainty);
P_before = Pmech_reduced;
%% ========================================================================
% 7.1 SEMPLIFICAZIONE DELLA RAPPRESENTAZIONE LFT
% ========================================================================
%
% Riduciamo le occorrenze ripetute delle incertezze mantenendo:
%
%   J_alpha
%   l
%
% come parametri incerti.
%
% La semplificazione non modifica nominalmente il modello e viene
% utilizzata solo per ottenere una rappresentazione LFT piu' compatta
% per la mu-synthesis.

Pmech_reduced = ...
    simplify( ...
        Pmech_reduced, ...
        'full');

fprintf('\n============================================================\n');
fprintf('INCERTEZZE DOPO SEMPLIFICAZIONE LFT\n');
fprintf('============================================================\n');

[~,~,blk_mech] = lftdata(Pmech_reduced);

for k = 1:numel(blk_mech)

    fprintf('%-12s -> %2d occorrenze\n', ...
        blk_mech(k).Name, ...
        blk_mech(k).Occurrences);

end

P_after = Pmech_reduced;
D = simplify(P_after - P_before,'full');

norm(D.NominalValue,inf)
[M1,Delta1] = lftdata(P_before);
[M2,Delta2] = lftdata(P_after);
fieldnames(P_before.Uncertainty)
fieldnames(P_after.Uncertainty)

fprintf('Incertezze PRIMA:\n');
disp(P_before.Uncertainty)

fprintf('Incertezze DOPO:\n');
disp(P_after.Uncertainty)

fprintf('Norma differenza nominale:\n');
Dnom = minreal( ...
    P_before.NominalValue - P_after.NominalValue, ...
    1e-8);

fprintf('||Pnom_before-Pnom_after||inf = %.6e\n', ...
    norm(Dnom,inf));

Ns = 500;

[P_before_s,SampleValues] = usample(P_before,Ns);

P_after_s = usubs(P_after,SampleValues);

err = zeros(Ns,1);

for k = 1:Ns
    err(k) = norm( ...
        P_before_s(:,:,k) - P_after_s(:,:,k), ...
        inf);
end

fprintf('Max errore sui %d campioni = %.6e\n', ...
    Ns,max(err));
%% ========================================================================
% 8. ESTRAZIONE DELLE USCITE ANGOLARI
% ========================================================================

C_angles = [
    1 0 0 0;
    0 0 1 0
];

Pq_mech_reduced = uss( ...
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
fprintf('\n============================================================\n');
fprintf('INCERTEZZE FINALI G UNCERTAIN LUMPED\n');
fprintf('============================================================\n');

disp(G_uncertain_lumped.Uncertainty);
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

G_uncertain_lumped_scaled = ...
    simplify( ...
        G_uncertain_lumped_scaled, ...
        'full');

G_uncertain_lumped_scaled.InputName = {
    'ubar1'
    'ubar2'
};

G_uncertain_lumped_scaled.OutputName = {
    'alpha_bar'
    'beta_bar'
};

fprintf('\n============================================================\n');
fprintf('LFT FINALE MODELLO PER MU-SYNTHESIS\n');
fprintf('============================================================\n');

[~,~,blk_final] = lftdata(G_uncertain_lumped_scaled);

for k = 1:numel(blk_final)

    fprintf('%-12s -> %2d occorrenze\n', ...
        blk_final(k).Name, ...
        blk_final(k).Occurrences);

end

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