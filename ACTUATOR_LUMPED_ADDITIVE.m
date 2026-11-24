%% ========================================================================
% MODELLO LUMPED ADDITIVO DEGLI ATTUATORI - ELICOTTERO 2-DOF
%
% SCOPO
%   Costruisce un secondo modello aggregato, additivo, senza sovrascrivere
%   ACTUATOR_LUMPED.mat (modello moltiplicativo già esistente).
%
% MODELLO
%   G_act_i = G_act_nom + WA_i * Delta_add_i
%
%   WA viene adattato all'inviluppo dell'ERRORE ASSOLUTO dell'intera
%   famiglia non-lumped degli attuatori, comprendendo:
%     - omega_n in +/-10%;
%     - tau_d in +/-20%;
%     - residuo W_Pade*Delta_Pade, con ||Delta_Pade||inf <= 1.
%
% PREREQUISITI
%   Eseguire prima dataset_elicottero.m e ACTUATOR_LUMPED.m nella stessa
%   sessione MATLAB, così G_act_nom, W_Pade, omega_n, zeta, tau_d_nom sono
%   disponibili. Il modello meccanico ridotto viene letto dal vecchio MAT.
%
% OUTPUT SEPARATO
%   ACTUATOR_LUMPED_ADDITIVE.mat
%
% Il vecchio ACTUATOR_LUMPED.mat e il vecchio MU_controller.mat NON vengono
% sovrascritti. Si creano anche copie di sicurezza se non già presenti.
% ========================================================================

close all;
clc;

fprintf('============================================================\n');
fprintf(' MODELLO LUMPED ADDITIVO DEGLI ATTUATORI\n');
fprintf('============================================================\n');

%% 0. CONTROLLI E COPIE DI SICUREZZA

if ~isfile('HINF_workspace.mat')
    error('HINF_workspace.mat non trovato. Eseguire prima HINF_SETUP.');
end
if ~isfile('ACTUATOR_LUMPED.mat')
    error('ACTUATOR_LUMPED.mat non trovato. Eseguire prima il setup lumped moltiplicativo.');
end

requiredVars = {'G_act_nom','W_Pade','omega_n','zeta','tau_d_nom'};
missingVars = {};
for k = 1:numel(requiredVars)
    if ~exist(requiredVars{k},'var')
        missingVars{end+1} = requiredVars{k}; %#ok<AGROW>
    end
end
if ~isempty(missingVars)
    error(['Variabili non presenti nel workspace: %s. Eseguire prima ' ...
           'dataset_elicottero.m e ACTUATOR_LUMPED.m nella stessa sessione.'], ...
           strjoin(missingVars, ', '));
end

% Copie di sicurezza del modello e del controllore moltiplicativi originali.
% Non sovrascrive le copie di sicurezza se esistono già.
if isfile('ACTUATOR_LUMPED.mat') && ~isfile('ACTUATOR_LUMPED_MULTIPLICATIVE_BACKUP.mat')
    copyfile('ACTUATOR_LUMPED.mat','ACTUATOR_LUMPED_MULTIPLICATIVE_BACKUP.mat');
    fprintf('Backup creato: ACTUATOR_LUMPED_MULTIPLICATIVE_BACKUP.mat\n');
end
if isfile('MU_controller.mat') && ~isfile('MU_controller_MULTIPLICATIVE_BACKUP.mat')
    copyfile('MU_controller.mat','MU_controller_MULTIPLICATIVE_BACKUP.mat');
    fprintf('Backup creato: MU_controller_MULTIPLICATIVE_BACKUP.mat\n');
end

% Importa soltanto le normalizzazioni e il modello meccanico già validato.
Ssetup = load('HINF_workspace.mat','Dy_inv','Du');
Dy_inv = Ssetup.Dy_inv;
Du = Ssetup.Du;
Smech = load('ACTUATOR_LUMPED.mat','Pq_mech_reduced','G_uncertain_lumped_scaled');
Pq_mech_reduced = Smech.Pq_mech_reduced;
G_uncertain_lumped_scaled_old = Smech.G_uncertain_lumped_scaled;

%% 1. INVILUPPO DELL'ERRORE ASSOLUTO DELL'ATTUATORE

% Parametri della famiglia coerenti con dataset_elicottero.m.
% La griglia copre omega_n +/-10% e tau_d +/-20%.
omega_fit = logspace(-2,4,2000);
jw_fit = 1i*omega_fit(:);
wn_fit = linspace(0.90*omega_n,1.10*omega_n,17);
tau_fit = linspace(0.80*tau_d_nom,1.20*tau_d_nom,41);

Gnom_resp = squeeze(freqresp(G_act_nom,omega_fit));
Gnom_resp = Gnom_resp(:);
Wp_abs_fit = squeeze(abs(freqresp(W_Pade,omega_fit)));
Wp_abs_fit = Wp_abs_fit(:);

E_add_env = zeros(numel(omega_fit),1);
for iwn = 1:numel(wn_fit)
    wn = wn_fit(iwn);
    G2nd_resp = wn^2 ./ ...
        (jw_fit.^2 + 2*zeta*wn*jw_fit + wn^2);

    for itau = 1:numel(tau_fit)
        tau = tau_fit(itau);
        Gpade_resp = (1 - tau*jw_fit/2) ./ (1 + tau*jw_fit/2);

        % Bound dell'errore assoluto completo:
        % |G2nd(wn)*Gpade(tau) - G_act_nom| + |G2nd(wn)|*|W_Pade|
        % Il secondo termine copre ogni Delta_Pade con norma Hinf <= 1.
        E_candidate = abs(G2nd_resp.*Gpade_resp - Gnom_resp) + ...
                      abs(G2nd_resp).*Wp_abs_fit;
        E_add_env = max(E_add_env,E_candidate);
    end
end

fprintf('\n============================================================\n');
fprintf(' INVILUPPO ERRORE ASSOLUTO ATTUATORI\n');
fprintf('============================================================\n');
fprintf('Combinazioni omega_n/tau_d analizzate: %d\n',numel(wn_fit)*numel(tau_fit));
fprintf('Massimo inviluppo errore assoluto = %.6e\n',max(E_add_env));

%% 2. FIT DI UN PESO ADDITIVO STABILE, PROPRIO E MINIMUM-PHASE

Margin_WA = 1.15;
E_add_fit = max(Margin_WA*E_add_env,1e-8);
E_add_frd = frd(reshape(E_add_fit,1,1,[]),omega_fit);

OrderWA = 6;
RelativeDegreeWA = 2;
CfitWA = struct();
CfitWA.LowerBound = E_add_frd;
CfitWA.UpperBound = Inf;

% Vincoliamo la copertura del modulo e imponiamo RD=2, coerente con il
% decadimento ad alta frequenza dell'errore assoluto dell'attuatore.
WA_fit = fitmagfrd(E_add_frd,OrderWA,RelativeDegreeWA,[],CfitWA);
WA = minreal(tf(WA_fit),1e-8);

fprintf('\nVerifica peso additivo WA:\n');
fprintf('Ordine WA = %d\n',order(WA));
fprintf('WA stabile = %d\n',isstable(WA));
fprintf('WA proprio = %d\n',isproper(WA));
if ~isstable(WA) || ~isproper(WA)
    error('WA non è stabile e proprio. Non procedere alla sintesi.');
end

%% 3. VALIDAZIONE SU GRIGLIA DENSA

omega_val = logspace(-2,4,4000);
jw_val = 1i*omega_val(:);
wn_val = linspace(0.90*omega_n,1.10*omega_n,31);
tau_val = linspace(0.80*tau_d_nom,1.20*tau_d_nom,61);

Gnom_val = squeeze(freqresp(G_act_nom,omega_val));
Gnom_val = Gnom_val(:);
Wp_abs_val = squeeze(abs(freqresp(W_Pade,omega_val)));
Wp_abs_val = Wp_abs_val(:);
mag_WA_val = squeeze(abs(freqresp(WA,omega_val)));
mag_WA_val = mag_WA_val(:);

E_add_env_val = zeros(numel(omega_val),1);
for iwn = 1:numel(wn_val)
    wn = wn_val(iwn);
    G2nd_resp = wn^2 ./ ...
        (jw_val.^2 + 2*zeta*wn*jw_val + wn^2);

    for itau = 1:numel(tau_val)
        tau = tau_val(itau);
        Gpade_resp = (1 - tau*jw_val/2) ./ (1 + tau*jw_val/2);
        E_candidate = abs(G2nd_resp.*Gpade_resp - Gnom_val) + ...
                      abs(G2nd_resp).*Wp_abs_val;
        E_add_env_val = max(E_add_env_val,E_candidate);
    end
end

ratio_WA_val = E_add_env_val ./ max(mag_WA_val,realmin);
ratio_WA_val_margin = Margin_WA*E_add_env_val ./ max(mag_WA_val,realmin);

fprintf('\n============================================================\n');
fprintf(' VALIDAZIONE PESO ADDITIVO SU GRIGLIA DENSA\n');
fprintf('============================================================\n');
fprintf('Combinazioni parametriche di validazione: %d\n',numel(wn_val)*numel(tau_val));
fprintf('Massimo rapporto errore/|WA| = %.6f\n',max(ratio_WA_val));
fprintf('Massimo rapporto con margine 15%% = %.6f\n',max(ratio_WA_val_margin));
fprintf('Frequenze fuori copertura base = %d / %d\n', ...
    sum(ratio_WA_val > 1+1e-6),numel(omega_val));
fprintf('Frequenze che violano margine = %d / %d\n', ...
    sum(ratio_WA_val_margin > 1+1e-6),numel(omega_val));

figure('Name','Lumped additivo - copertura errore assoluto');
semilogx(omega_val,E_add_env_val,'LineWidth',1.5); hold on;
semilogx(omega_val,mag_WA_val,'LineWidth',1.5);
semilogx(omega_val,Margin_WA*E_add_env_val,'--','LineWidth',1.2);
grid on;
xlabel('\omega [rad/s]'); ylabel('Modulo');
title('Inviluppo dell''errore assoluto e peso additivo W_A');
legend('Inviluppo errore assoluto','|W_A|','1.15 x inviluppo','Location','best');

% Non lanciare una sintesi se la copertura numerica del peso additivo fallisce.
if max(ratio_WA_val_margin) > 1+1e-3
    error(['Il fit WA non copre l''inviluppo con il margine del 15%% sulla ' ...
           'griglia densa. Aumentare OrderWA (per esempio a 6) e rieseguire.']);
end

%% 4. MODELLO ADDITIVO LUMPED, DUE ATTUATORI INDIPENDENTI

WA1 = WA;
WA2 = WA;
Delta_add1 = ultidyn('Delta_add1',[1 1]);
Delta_add2 = ultidyn('Delta_add2',[1 1]);

G1_unc_add_lumped = G_act_nom + WA1*Delta_add1;
G2_unc_add_lumped = G_act_nom + WA2*Delta_add2;

Gact_additive_lumped = blkdiag(G1_unc_add_lumped,G2_unc_add_lumped);
Gact_additive_lumped.InputName = {'u1';'u2'};
Gact_additive_lumped.OutputName = {'F1';'F2'};

Pq_mech_reduced.InputName = {'F1';'F2'};
Pq_mech_reduced.OutputName = {'alpha';'beta'};
G_uncertain_additive_lumped = simplify(Pq_mech_reduced*Gact_additive_lumped,'full');
G_uncertain_additive_lumped.InputName = {'u1';'u2'};
G_uncertain_additive_lumped.OutputName = {'alpha';'beta'};

G_uncertain_additive_lumped_scaled = simplify( ...
    Dy_inv*G_uncertain_additive_lumped*Du,'full');
G_uncertain_additive_lumped_scaled.InputName = {'ubar1';'ubar2'};
G_uncertain_additive_lumped_scaled.OutputName = {'alpha_bar';'beta_bar'};

% Controllo che il nuovo modello abbia lo stesso nominale del vecchio lumped.
Gnom_add = G_uncertain_additive_lumped_scaled.NominalValue;
Gnom_old = G_uncertain_lumped_scaled_old.NominalValue;
nominalDifference = minreal(Gnom_add-Gnom_old,1e-8);
nominalDifferenceNorm = norm(nominalDifference,inf);
fprintf('\nDifferenza fra nominale lumped moltiplicativo e additivo = %.6e\n', ...
    nominalDifferenceNorm);

fprintf('\nBlocchi incerti del nuovo modello lumped additivo:\n');
[~,~,blk_add] = lftdata(G_uncertain_additive_lumped_scaled);
for k = 1:numel(blk_add)
    fprintf('%-14s -> %2d occorrenze\n',blk_add(k).Name,blk_add(k).Occurrences);
end

%% 5. SALVATAGGIO SEPARATO - NON SOVRASCRIVE IL VECCHIO FILE

outFile = 'ACTUATOR_LUMPED_ADDITIVE.mat';
if isfile(outFile)
    backupName = ['ACTUATOR_LUMPED_ADDITIVE_previous_' ...
        char(datetime('now','Format','yyyyMMdd_HHmmss')) '.mat'];
    copyfile(outFile,backupName);
    fprintf('Versione additiva precedente salvata in %s\n',backupName);
end

save(outFile, ...
    'WA','WA1','WA2','Delta_add1','Delta_add2', ...
    'Gact_additive_lumped','G_uncertain_additive_lumped', ...
    'G_uncertain_additive_lumped_scaled','Pq_mech_reduced', ...
    'E_add_env','E_add_env_val','omega_fit','omega_val', ...
    'wn_fit','tau_fit','wn_val','tau_val', ...
    'ratio_WA_val','ratio_WA_val_margin', ...
    'Margin_WA','OrderWA','RelativeDegreeWA', ...
    'G_act_nom','W_Pade','omega_n','zeta','tau_d_nom', ...
    'nominalDifferenceNorm');

fprintf('\n============================================================\n');
fprintf('MODELLO ADDITIVO SALVATO IN %s\n',outFile);
fprintf('Il file ACTUATOR_LUMPED.mat non è stato modificato.\n');
fprintf('============================================================\n');

%% TEST DI INCLUSIONE NEL MODELLO ADDITIVO

Sadd = load('ACTUATOR_LUMPED_ADDITIVE.mat','WA');
WA = Sadd.WA;

s = tf('s');

wn_test = omega_n;
tau_test = 1.20*tau_d_nom;

G2_test = wn_test^2 / ...
    (s^2 + 2*zeta*wn_test*s + wn_test^2);

Gpade_test = ...
    (1 - tau_test*s/2) / ...
    (1 + tau_test*s/2);

% Residuo Delta_Pade posto a zero
G_NL_test = G2_test*Gpade_test;

% Perturbazione additiva necessaria
Delta_add_required = minreal( ...
    (G_NL_test - G_act_nom)/WA, 1e-8);

fprintf('\nTEST INCLUSIONE LUMPED ADDITIVO\n');
fprintf('Delta richiesta stabile = %d\n', ...
    isstable(Delta_add_required));
fprintf('Delta richiesta propria = %d\n', ...
    isproper(Delta_add_required));

if isstable(Delta_add_required) && isproper(Delta_add_required)
    gamma_add_required = hinfnorm(Delta_add_required);
    fprintf('||Delta richiesta||inf = %.6f\n', gamma_add_required);
    fprintf('Bound unitario rispettato = %d\n', ...
        gamma_add_required <= 1+1e-6);
end

%% VERIFICA DEL PESO ADDITIVO E DEI CASI PARAMETRICI ESTREMI

S = load('ACTUATOR_LUMPED_ADDITIVE.mat', ...
    'WA','G_act_nom','omega_n','zeta','tau_d_nom');

WA = S.WA;
G_act_nom = S.G_act_nom;

zWA = zero(WA);

fprintf('\n============================================\n');
fprintf('VERIFICA STRUTTURALE DEL MODELLO ADDITIVO\n');
fprintf('============================================\n');

fprintf('WA stabile = %d\n',isstable(WA));
fprintf('WA proprio = %d\n',isproper(WA));
fprintf('WA minimum-phase = %d\n',all(real(zWA)<-1e-8));

disp('Zeri di WA:');
disp(zWA);

% Verifica i quattro angoli della famiglia parametrica:
% wn = +/-10%, tau_d = +/-20%.
wn_vec = [0.90 1.10]*S.omega_n;
tau_vec = [0.80 1.20]*S.tau_d_nom;

fprintf('\nCasi parametrici estremi:\n');

for i = 1:numel(wn_vec)
    for j = 1:numel(tau_vec)

        wn = wn_vec(i);
        tau = tau_vec(j);
        s = tf('s');

        G2 = wn^2/(s^2 + 2*S.zeta*wn*s + wn^2);
        Gpade = (1-tau*s/2)/(1+tau*s/2);
        G_NL_test = G2*Gpade;

        Delta_req = minreal( ...
            (G_NL_test-G_act_nom)/WA,1e-8);

        if isstable(Delta_req) && isproper(Delta_req)
            gamma_req = hinfnorm(Delta_req);
        else
            gamma_req = Inf;
        end

        fprintf(['wn = %6.2f | tau = %.4f | stabile = %d | ' ...
            'propria = %d | norma = %.6f\n'], ...
            wn,tau,isstable(Delta_req), ...
            isproper(Delta_req),gamma_req);
    end
end