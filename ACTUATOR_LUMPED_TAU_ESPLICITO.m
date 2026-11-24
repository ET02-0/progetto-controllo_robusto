%% ================================================================
% LUMPED MOLTIPLICATIVO CON TAU_D ESPLICITO
%
% Non sovrascrive ACTUATOR_LUMPED.mat o MU_controller.mat.
%
% Espliciti:
%   tau_d1, tau_d2: +/-20%
%
% Lumped moltiplicativamente:
%   omega_n +/-10% e residuo dell'approssimazione di Padé.
% ================================================================

close all;
clc;

%% 0. CONTROLLI E CARICAMENTO

requiredVars = {'G_act_nom','W_Pade','omega_n','zeta','tau_d_nom'};

for k = 1:numel(requiredVars)
    if ~exist(requiredVars{k},'var')
        error(['Variabile %s mancante. Eseguire prima ' ...
               'dataset_elicottero.m.'], requiredVars{k});
    end
end

assert(isfile('HINF_workspace.mat'), ...
    'HINF_workspace.mat non trovato.');

assert(isfile('ACTUATOR_LUMPED.mat'), ...
    'ACTUATOR_LUMPED.mat non trovato.');

H = load('HINF_workspace.mat','Du','Dy_inv');

M = load('ACTUATOR_LUMPED.mat', ...
    'Pq_mech_reduced', ...
    'G_uncertain_lumped_scaled');

Du = H.Du;
Dy_inv = H.Dy_inv;
Pq_mech_reduced = M.Pq_mech_reduced;

%% 1. CALCOLO DEL PESO MOLTIPLICATIVO RESIDUO
%
% Il ritardo tau_d NON viene incluso nell'inviluppo da approssimare:
% resta esplicitamente nel modello nominale parametrico.
%
% Il peso residuo copre:
%   - variazione di omega_n;
%   - residuo di Padé già descritto da W_Pade.
%
% Rapporto rispetto al secondo ordine nominale:
% Rwn = G2(omega_n) / G2_nom
%
% Bound moltiplicativo:
% |Rwn - 1| + |Rwn|*|W_Pade|

omega = logspace(-2,4,2000);
jw = 1i*omega(:);

wn_vec = linspace(0.90*omega_n,1.10*omega_n,17);

G2_nom = tf(omega_n^2, ...
    [1, 2*zeta*omega_n, omega_n^2]);

G2_nom_resp = squeeze(freqresp(G2_nom,omega));
G2_nom_resp = G2_nom_resp(:);

Wp_abs = squeeze(abs(freqresp(W_Pade,omega)));
Wp_abs = Wp_abs(:);

E_env = zeros(numel(omega),1);

for iwn = 1:numel(wn_vec)

    wn = wn_vec(iwn);

    G2_wn_resp = wn^2 ./ ...
        (jw.^2 + 2*zeta*wn*jw + wn^2);

    Rwn = G2_wn_resp ./ G2_nom_resp;

    E_candidate = abs(Rwn - 1) + abs(Rwn).*Wp_abs;

    E_env = max(E_env,E_candidate);

end

fprintf('\n============================================\n');
fprintf('INVILUPPO MOLTIPLICATIVO RESIDUO\n');
fprintf('============================================\n');
fprintf('Massimo inviluppo = %.6f\n',max(E_env));

%% 2. FIT DEL PESO

Margin_WI_rest = 1.15;
OrderWt_rest = 4;

E_fit = max(Margin_WI_rest*E_env,1e-8);

E_fit_frd = frd(reshape(E_fit,1,1,[]),omega);

Cfit_rest = struct();
Cfit_rest.LowerBound = E_fit_frd;
Cfit_rest.UpperBound = Inf;

WI_rest_fit = fitmagfrd( ...
    E_fit_frd,OrderWt_rest,[],[],Cfit_rest);

WI_rest = minreal(tf(WI_rest_fit),1e-7);

fprintf('\nVerifica WI_rest:\n');
fprintf('Ordine = %d\n',order(WI_rest));
fprintf('Stabile = %d\n',isstable(WI_rest));
fprintf('Proprio = %d\n',isproper(WI_rest));

if ~isstable(WI_rest) || ~isproper(WI_rest)
    error('WI_rest non stabile o non proprio.');
end

%% 3. VALIDAZIONE DEL PESO SU GRIGLIA DENSA

omega_val = logspace(-2,4,4000);
jw_val = 1i*omega_val(:);

wn_val = linspace(0.90*omega_n,1.10*omega_n,31);

G2_nom_val = squeeze(freqresp(G2_nom,omega_val));
G2_nom_val = G2_nom_val(:);

Wp_abs_val = squeeze(abs(freqresp(W_Pade,omega_val)));
Wp_abs_val = Wp_abs_val(:);

WI_abs_val = squeeze(abs(freqresp(WI_rest,omega_val)));
WI_abs_val = WI_abs_val(:);

E_env_val = zeros(numel(omega_val),1);

for iwn = 1:numel(wn_val)

    wn = wn_val(iwn);

    G2_wn_val = wn^2 ./ ...
        (jw_val.^2 + 2*zeta*wn*jw_val + wn^2);

    Rwn_val = G2_wn_val ./ G2_nom_val;

    E_candidate_val = ...
        abs(Rwn_val - 1) + abs(Rwn_val).*Wp_abs_val;

    E_env_val = max(E_env_val,E_candidate_val);

end

ratio_val = E_env_val ./ max(WI_abs_val,realmin);
ratio_val_margin = ...
    Margin_WI_rest*E_env_val ./ max(WI_abs_val,realmin);

fprintf('\n============================================\n');
fprintf('VALIDAZIONE WI_rest\n');
fprintf('============================================\n');
fprintf('Massimo rapporto errore/WI = %.6f\n',max(ratio_val));
fprintf('Massimo rapporto con margine = %.6f\n', ...
    max(ratio_val_margin));
fprintf('Frequenze fuori copertura = %d / %d\n', ...
    sum(ratio_val > 1+1e-6),numel(omega_val));

figure('Name','Peso moltiplicativo residuo');
semilogx(omega_val,E_env_val,'LineWidth',1.5);
hold on;
semilogx(omega_val,WI_abs_val,'LineWidth',1.5);
grid on;
xlabel('\omega [rad/s]');
ylabel('Modulo');
legend('Inviluppo residuo','|WI_{rest}|','Location','best');
title('Copertura del peso moltiplicativo residuo');

if max(ratio_val_margin) > 1+1e-3
    error(['WI_rest non copre l''inviluppo con il margine richiesto. ' ...
        'Non procedere alla sintesi finché il fit non viene corretto.']);
end

%% 4. RITARDI ESPLICITI E INDIPENDENTI

s_tau = tf('s');

tau_d1 = ureal('tau_d1',tau_d_nom,'Percentage',20);
tau_d2 = ureal('tau_d2',tau_d_nom,'Percentage',20);

Gdelay_tau1 = ...
    (1 - tau_d1*s_tau/2)/(1 + tau_d1*s_tau/2);

Gdelay_tau2 = ...
    (1 - tau_d2*s_tau/2)/(1 + tau_d2*s_tau/2);

%% 5. INCERTEZZA MOLTIPLICATIVA RESIDUA

Delta_act1_tau = ultidyn('Delta_act1_tau',[1 1]);
Delta_act2_tau = ultidyn('Delta_act2_tau',[1 1]);

G1_unc_tau = ...
    G2_nom*Gdelay_tau1*(1 + WI_rest*Delta_act1_tau);

G2_unc_tau = ...
    G2_nom*Gdelay_tau2*(1 + WI_rest*Delta_act2_tau);

%% ================================================================
% SEMPLIFICAZIONE OCCORRENZE E CONTROLLO LFT
% ================================================================

%% 1. SEMPLIFICAZIONE DEI DUE ATTUATORI

Gact_lumped_tau_explicit = blkdiag(G1_unc_tau,G2_unc_tau);

Gact_lumped_tau_explicit = simplify( ...
    Gact_lumped_tau_explicit,'full');

Gact_lumped_tau_explicit.InputName = {'u1';'u2'};
Gact_lumped_tau_explicit.OutputName = {'F1';'F2'};

fprintf('\n============================================\n');
fprintf('LFT ATTUATORI DOPO SEMPLIFICAZIONE\n');
fprintf('============================================\n');

[~,~,blk_act_tau] = lftdata(Gact_lumped_tau_explicit);

for k = 1:numel(blk_act_tau)
    fprintf('%-24s -> %2d occorrenze\n', ...
        blk_act_tau(k).Name,blk_act_tau(k).Occurrences);
end

%% 2. CONNESSIONE CON LA MECCANICA GIÀ LUMPED

G_uncertain_lumped_tau_explicit = ...
    Pq_mech_reduced * Gact_lumped_tau_explicit;

G_uncertain_lumped_tau_explicit = simplify( ...
    G_uncertain_lumped_tau_explicit,'full');

G_uncertain_lumped_tau_explicit.InputName = {'u1';'u2'};
G_uncertain_lumped_tau_explicit.OutputName = {'alpha';'beta'};

%% 3. NORMALIZZAZIONE E SEMPLIFICAZIONE FINALE

G_uncertain_lumped_tau_explicit_scaled = simplify( ...
    Dy_inv * G_uncertain_lumped_tau_explicit * Du,'full');

G_uncertain_lumped_tau_explicit_scaled.InputName = {'ubar1';'ubar2'};
G_uncertain_lumped_tau_explicit_scaled.OutputName = ...
    {'alpha_bar';'beta_bar'};

fprintf('\n============================================\n');
fprintf('LFT FINALE DOPO SEMPLIFICAZIONE\n');
fprintf('============================================\n');

[~,~,blk_tau_final] = ...
    lftdata(G_uncertain_lumped_tau_explicit_scaled);

for k = 1:numel(blk_tau_final)
    fprintf('%-24s -> %2d occorrenze\n', ...
        blk_tau_final(k).Name,blk_tau_final(k).Occurrences);
end

%% 4. VERIFICA PRESENZA DEI DUE RITARDI ESPLICITI

uncNames = fieldnames( ...
    G_uncertain_lumped_tau_explicit_scaled.Uncertainty);

disp('Nomi dei blocchi incerti presenti nel modello finale:');
disp(uncNames);

if ~any(strcmp(uncNames,'tau_d1')) || ...
        ~any(strcmp(uncNames,'tau_d2'))
    warning(['Controllare la rappresentazione LFT: ' ...
        'tau_d1 e tau_d2 non compaiono entrambi nei blocchi incerti.']);
end
%% 9. SALVATAGGIO SEPARATO

save('ACTUATOR_LUMPED_TAU_ESPLICITO.mat', ...
    'WI_rest','Delta_act1_tau','Delta_act2_tau', ...
    'tau_d1','tau_d2', ...
    'Gact_lumped_tau_explicit', ...
    'G_uncertain_lumped_tau_explicit', ...
    'G_uncertain_lumped_tau_explicit_scaled', ...
    'Pq_mech_reduced', ...
    'E_env','E_env_val','omega','omega_val', ...
    'ratio_val','ratio_val_margin', ...
    'Margin_WI_rest','OrderWt_rest');

fprintf('\nSalvato ACTUATOR_LUMPED_TAU_ESPLICITO.mat\n');
fprintf('I file precedenti non sono stati sovrascritti.\n');