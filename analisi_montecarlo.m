%% ==========================================
% MONTE CARLO ANALITICO (SENZA SIMULINK)
% Elicottero 2-DOF - Validazione Robusta
% ==========================================
close all; clc;

disp('--- Avvio Analisi Monte Carlo Analitica ---')

%% 1. CARICAMENTO DATI E CONTROLLORI

load('HINF_workspace.mat');

load('HINF_controllers.mat', ...
    'K_mix','K_hinfsyn','K_pidcomp');

% Carica il controllore MU additivo selezionato e validato
fileMUadd = 'MU_controller_ADDITIVE_15PCT_COLDSTART_MINREAL.mat';

assert(isfile(fileMUadd), ...
    'File del controllore MU additivo non trovato.');

Madd = load(fileMUadd, ...
    'K_mu_additive_phys', ...
    'controllerSelected', ...
    'wc_add_on_NL', ...
    'stabmarg_add_NL');

assert(isfield(Madd,'K_mu_additive_phys') && ...
    ~isempty(Madd.K_mu_additive_phys), ...
    'K_mu_additive_phys non disponibile.');

K_mu = ss(Madd.K_mu_additive_phys);

fprintf('Controllore MU selezionato: %s\n', ...
    Madd.controllerSelected);

fprintf('Ordine controllore MU additivo = %d\n',order(K_mu));

fprintf('RP upper sul modello non-lumped = %.6f\n', ...
    Madd.wc_add_on_NL.UpperBound);

fprintf('RS lower sul modello non-lumped = %.6f\n', ...
    Madd.stabmarg_add_NL.LowerBound);
%% VERIFICA MODELLO INCERTO COMPLETO

% Incertezze che devono essere presenti nel modello utilizzato
unc_richieste = { ...
    'J_alpha', 'Jy', 'Jz', 'l', 'm', ...
    'eps_p', 'eps_y', ...
    'omega_n1', 'omega_n2', ...
    'tau_d1', 'tau_d2', ...
    'Delta_Pade1', 'Delta_Pade2'};

% Nomi delle incertezze effettivamente presenti
unc_presenti = fieldnames(G_uncertain.Uncertainty);

% Verifica delle incertezze richieste
unc_mancanti = setdiff(unc_richieste, unc_presenti);

if ~isempty(unc_mancanti)
    error('Mancano questi blocchi di incertezza: %s', ...
        strjoin(unc_mancanti, ', '));
end

fprintf('Verifica delle incertezze completata: tutti i blocchi richiesti sono presenti.\n');

% Scegli il controllore da testare (es. K_mix, K_hinfsyn, K_pidcomp o K_mu)
controllers = {
    K_mix
    K_hinfsyn
    K_pidcomp
    K_mu
    };

controllerNames = {
    'mixsyn'
    'hinfsyn'
    'PID+comp'
    'mu-synthesis additiva'
    };

Nc = numel(controllers);

N_campioni = 500;
rng('default');

%% 2. COSTRUZIONE MODELLO DISTURBI (Gd_uncertain)
% Il setup H-inf base contiene G_uncertain (da u a y).
% Dobbiamo creare Gd_uncertain (da d a y) per testare i disturbi aerodinamici.
% Usiamo le matrici A_unc e Bd_unc generate dal dataset_elicottero.m
if ~exist('A_unc', 'var')
    error('Esegui prima dataset_elicottero.m per avere le matrici di stato in workspace.');
end

C_angles = [1 0 0 0; 0 0 1 0]; % Estrazione di [alpha; beta]
D_d = zeros(2,2);
Gd_uncertain = uss(A_unc, Bd_unc, C_angles, D_d);

% Uniamo impianto e disturbi in un unico oggetto per campionarli insieme in modo coerente
Gall_unc = [G_uncertain, Gd_uncertain];
Gall_samples = usample(Gall_unc, N_campioni);

%% 3. PARAMETRI DI TEST (Dal dataset)
amp_step_alpha = deg2rad(1);
amp_step_beta  = deg2rad(2);
amp_dist_alpha = 5e-3;       % Ampiezza disturbo pitch [N*m]
amp_dist_beta  = 2e-3;       % Ampiezza disturbo yaw [N*m]

% Limite di saturazione attuatori (variazioni rispetto all'equilibrio)
umax = 2.5; 

% Soglia e orizzonte per il calcolo del tempo di assestamento
settlingThreshold = 0.05;
% Vettore dei tempi per la simulazione analitica
t_sim = 0:0.005:30; 

%% PARAMETRI PER LA VERIFICA DI RECUPERO DAI DISTURBI
t_disturbo = 10;                  % Istante di applicazione [s]
banda_rec_alpha = deg2rad(0.3);   % Banda di errore pitch [rad]
banda_rec_beta  = deg2rad(0.7);   % Banda di errore yaw [rad]
durata_banda = 1;                 % Permanenza richiesta nella banda [s]
% Soglie RMS negli ultimi 5 secondi
limite_RMS_alpha = deg2rad(0.3);
limite_RMS_beta  = deg2rad(0.7);

dt_sim = t_sim(2) - t_sim(1);

idx_disturbo = find(t_sim >= t_disturbo, 1, 'first');
assert(~isempty(idx_disturbo) && ...
    abs(t_sim(idx_disturbo) - t_disturbo) <= dt_sim/2, ...
    'L''istante del disturbo non coincide con un campione di t_sim.');

idx_ultimi5 = t_sim >= t_sim(end) - 5;

%% 4. INIZIALIZZAZIONE METRICHE

Ntot = Nc * N_campioni;

resultController = strings(Ntot,1);
resultSample     = zeros(Ntot,1);

Stable = false(Ntot,1);

PitchSettling = nan(Ntot,1);
YawSettling   = nan(Ntot,1);

PitchSettledWithinHorizon = false(Ntot,1);
YawSettledWithinHorizon   = false(Ntot,1);

PitchBW = nan(Ntot,1);
YawBW   = nan(Ntot,1);

PitchOvershoot = nan(Ntot,1);
YawOvershoot   = nan(Ntot,1);

PitchSSerror = nan(Ntot,1);
YawSSerror   = nan(Ntot,1);

YawFromPitchPeak = nan(Ntot,1);
PitchFromYawPeak = nan(Ntot,1);

PitchDistPeak = nan(Ntot,1);
YawDistPeak   = nan(Ntot,1);


max_u1_cmd = nan(Ntot,1);
max_u2_cmd = nan(Ntot,1);

min_u1_cmd = nan(Ntot,1);
max_u1_signed = nan(Ntot,1);

min_u2_cmd = nan(Ntot,1);
max_u2_signed = nan(Ntot,1);

%% METRICHE DI RECUPERO DAI DISTURBI

PitchRecoveryTime = nan(Ntot,1);
YawRecoveryTime   = nan(Ntot,1);

% true se l'errore resta nella banda per almeno 1 secondo
PitchRecovered1s = false(Ntot,1);
YawRecovered1s   = false(Ntot,1);

% Rispetto dei limiti temporali di recupero
PitchRecoveryPass = false(Ntot,1);
YawRecoveryPass   = false(Ntot,1);

% Errore residuo massimo e RMS negli ultimi 5 secondi
PitchResidualPeakLast5s = nan(Ntot,1);
YawResidualPeakLast5s   = nan(Ntot,1);

PitchResidualRMSLast5s = nan(Ntot,1);
YawResidualRMSLast5s   = nan(Ntot,1);

% Verifica conformità RMS negli ultimi 5 secondi
PitchResidualRMSPass = false(Ntot,1);
YawResidualRMSPass   = false(Ntot,1);

%% LIMITI FISICI DEGLI ATTUATORI (N)
u1_min_lim = -2.5;
u1_max_lim =  2.5;

u2_min_lim = -2.5;
u2_max_lim =  2.5;

row = 0;
%% 5. CICLO MONTE CARLO
for ic = 1:Nc
    K_test = controllers{ic};
    for i = 1:N_campioni
        row = row + 1;

        resultController(row) = string(controllerNames{ic});
        resultSample(row) = i;
        
        % Estrazione del singolo campione
        G_sample = Gall_samples(:,:,i);
        Gu = G_sample(:, 1:2); % Sottomatrice da u a y (Insegumento)
        Gd = G_sample(:, 3:4); % Sottomatrice da d a y (Disturbi)
        
        % Costruzione delle funzioni di anello chiuso
        L  = Gu * K_test;
        S  = feedback(eye(2), L);
        T  = feedback(L, eye(2));
        KS = K_test * S; % Funzione di sensitività del controllo (da r a u)
        Td = S * Gd;     % Funzione di sensitività ai disturbi (da d a y)
        
        % --- STABILITÀ ---
        Stable(row) = all(real(pole(T)) < 0);
        if ~Stable(row)
            fprintf('%s - Campione %d: INSTABILE!\n', ...
                controllerNames{ic}, i);
            continue; 
        end
        
        % --- BANDA PASSANTE ---
        PitchBW(row) = bandwidth(T(1,1));
        YawBW(row)   = bandwidth(T(2,2));
        
        %% ---------------------------------------------------------------
        % TRACKING PITCH
        % ---------------------------------------------------------------
        
        [y_track, ~] = step(T(:,1) * amp_step_alpha, t_sim);
        y_track = squeeze(y_track);
        
        [u_track, ~] = step(KS(:,1) * amp_step_alpha, t_sim);
        u_track = squeeze(u_track);
        
        % Valore finale reale
        yss_alpha = dcgain(T(1,1)) * amp_step_alpha;
        
        % Tempo di assestamento al 5%
        infoPitch = stepinfo( ...
            y_track(:,1), t_sim, yss_alpha, ...
            'SettlingTimeThreshold', settlingThreshold);

        PitchSettling(row) = infoPitch.SettlingTime;
        PitchOvershoot(row) = infoPitch.Overshoot;

        if ~isfinite(PitchSettling(row))
            [PitchSettling(row), PitchSettledWithinHorizon(row)] = ...
                settlingTimeFromTrace( ...
                t_sim, y_track(:,1), yss_alpha, settlingThreshold);
        else
            PitchSettledWithinHorizon(row) = true;
        end

        % Errore a regime
        PitchSSerror(row) = ...
            100 * abs(amp_step_alpha - yss_alpha) / abs(amp_step_alpha);
        
        % Cross-coupling alpha -> beta
        YawFromPitchPeak(row) = max(abs(y_track(:,2)));
            
        
        %% SFORZO ATTUATORI SU ENTRAMBI I RIFERIMENTI

        % Risposta dei comandi al gradino di yaw
        [u_yaw,~] = step(KS(:,2)*amp_step_beta,t_sim);
        u_yaw = squeeze(u_yaw);

        % Comandi generati dai due gradini separati
        U1_all = [u_track(:,1), u_yaw(:,1)];
        U2_all = [u_track(:,2), u_yaw(:,2)];

        % Picchi assoluti
        max_u1_cmd(row) = max(abs(U1_all(:)));
        max_u2_cmd(row) = max(abs(U2_all(:)));

        % Estremi con segno, necessari per i limiti asimmetrici
        min_u1_cmd(row) = min(U1_all(:));
        max_u1_signed(row) = max(U1_all(:));

        min_u2_cmd(row) = min(U2_all(:));
        max_u2_signed(row) = max(U2_all(:));

        %% ---------------------------------------------------------------
        % TRACKING YAW
        % ---------------------------------------------------------------
        
        [y_yaw, ~] = step(T(:,2) * amp_step_beta, t_sim);
        y_yaw = squeeze(y_yaw);
        
        % Valore finale reale
        yss_beta = dcgain(T(2,2)) * amp_step_beta;
        
        % Tempo di assestamento al 5%
        infoYaw = stepinfo( ...
            y_yaw(:,2), t_sim, yss_beta, ...
            'SettlingTimeThreshold', settlingThreshold);

        YawSettling(row) = infoYaw.SettlingTime;
        YawOvershoot(row) = infoYaw.Overshoot;

        if ~isfinite(YawSettling(row))
            [YawSettling(row), YawSettledWithinHorizon(row)] = ...
                settlingTimeFromTrace( ...
                t_sim, y_yaw(:,2), yss_beta, settlingThreshold);
        else
            YawSettledWithinHorizon(row) = true;
        end
        
        % Errore a regime
        YawSSerror(row) = ...
            100 * abs(amp_step_beta - yss_beta) / abs(amp_step_beta);
        
        % Cross-coupling beta -> alpha
        PitchFromYawPeak(row) = max(abs(y_yaw(:,1)));
        
        %% ---------------------------------------------------------------
        % DISTURBI AERODINAMICI E RECUPERO DEL RIFERIMENTO
        % ---------------------------------------------------------------

        % Risposte ai gradini di disturbo, calcolate a partire da tau = 0.
        [y_da, ~] = step(Td(:,1) * amp_dist_alpha, t_sim);
        [y_db, ~] = step(Td(:,2) * amp_dist_beta,  t_sim);

        y_da = squeeze(y_da);
        y_db = squeeze(y_db);

        % Normalizza la disposizione delle matrici: tempo x 2 uscite
        if size(y_da,1) ~= numel(t_sim)
            y_da = y_da.';
        end
        if size(y_db,1) ~= numel(t_sim)
            y_db = y_db.';
        end

        % Superposizione dei due disturbi, applicati simultaneamente a t = 10 s.
        % Le risposte calcolate da step vengono traslate all'istante di disturbo.
        N_dopo = numel(t_sim) - idx_disturbo + 1;

        y_dist_tot = zeros(numel(t_sim), 2);

        y_dist_tot(idx_disturbo:end,:) = ...
            y_da(1:N_dopo,:) + y_db(1:N_dopo,:);

        % Risposta ai due riferimenti mantenuti durante il disturbo.
        % Il sistema è LTI, quindi le due risposte possono essere sommate.
        y_ref_tot = y_track + y_yaw;

        % Risposta complessiva: riferimenti + disturbi
        y_tot_dist = y_ref_tot + y_dist_tot;

        % Errori rispetto ai riferimenti richiesti
        e_alpha_dist = y_tot_dist(:,1) - amp_step_alpha;
        e_beta_dist  = y_tot_dist(:,2) - amp_step_beta;

        % Picco della sola risposta ai disturbi, dopo t = 10 s
        PitchDistPeak(row) = max( ...
            abs(y_dist_tot(idx_disturbo:end,1)));

        YawDistPeak(row) = max( ...
            abs(y_dist_tot(idx_disturbo:end,2)));

        % Tempo di recupero: primo rientro nella banda seguito da almeno
        % un secondo interamente trascorso all'interno della banda.
        [PitchRecoveryTime(row), PitchRecovered1s(row)] = ...
            recoveryTimeFromError( ...
            t_sim, e_alpha_dist, t_disturbo, ...
            banda_rec_alpha, durata_banda);

        [YawRecoveryTime(row), YawRecovered1s(row)] = ...
            recoveryTimeFromError( ...
            t_sim, e_beta_dist, t_disturbo, ...
            banda_rec_beta, durata_banda);

        % Verifica delle specifiche temporali di recupero
        PitchRecoveryPass(row) = ...
            PitchRecovered1s(row) && ...
            PitchRecoveryTime(row) <= 2;

        YawRecoveryPass(row) = ...
            YawRecovered1s(row) && ...
            YawRecoveryTime(row) <= 3;

        % Errore residuo negli ultimi 5 secondi della simulazione
        e_alpha_finale = e_alpha_dist(idx_ultimi5);
        e_beta_finale  = e_beta_dist(idx_ultimi5);

        PitchResidualPeakLast5s(row) = ...
            max(abs(e_alpha_finale));

        YawResidualPeakLast5s(row) = ...
            max(abs(e_beta_finale));

        PitchResidualRMSLast5s(row) = ...
            sqrt(mean(e_alpha_finale.^2));

        YawResidualRMSLast5s(row) = ...
            sqrt(mean(e_beta_finale.^2));

        % Verifica delle soglie RMS
        PitchResidualRMSPass(row) = ...
            PitchResidualRMSLast5s(row) <= limite_RMS_alpha;

        YawResidualRMSPass(row) = ...
            YawResidualRMSLast5s(row) <= limite_RMS_beta;
        % Stampa a schermo progresso
        fprintf('%s - Campione %d/%d elaborato.\n', ...
            controllerNames{ic}, i, N_campioni);
    end
end
%% VERIFICA SPECIFICHE SUL PICCO DELLA RISPOSTA AI DISTURBI

limite_picco_alpha = deg2rad(0.85);  % [rad]
limite_picco_beta  = deg2rad(5.0);  % [rad]

fprintf('\n============================================\n');
fprintf(' VERIFICA PICCO RISPOSTA AI DISTURBI\n');
fprintf('============================================\n');

for ic = 1:Nc

    idxCtrl = strcmp(resultController, ...
        string(controllerNames{ic})) & Stable;

    % Considera solo campioni stabili con picco disponibile
    validAlpha = idxCtrl & isfinite(PitchDistPeak);
    validBeta  = idxCtrl & isfinite(YawDistPeak);

    nAlpha = sum(validAlpha);
    nBeta  = sum(validBeta);

    % Campioni conformi e non conformi
    passAlpha = sum(validAlpha & ...
        PitchDistPeak <= limite_picco_alpha);
    failAlpha = sum(validAlpha & ...
        PitchDistPeak > limite_picco_alpha);

    passBeta = sum(validBeta & ...
        YawDistPeak <= limite_picco_beta);
    failBeta = sum(validBeta & ...
        YawDistPeak > limite_picco_beta);

    fprintf('\nControllore: %s\n', controllerNames{ic});

    fprintf(['Picco alpha: conformi %d/%d (%.1f%%), ' ...
        'superamenti %d/%d\n'], ...
        passAlpha, nAlpha, ...
        100*passAlpha/max(nAlpha,1), failAlpha, nAlpha);

    fprintf(['Picco beta:  conformi %d/%d (%.1f%%), ' ...
        'superamenti %d/%d\n'], ...
        passBeta, nBeta, ...
        100*passBeta/max(nBeta,1), failBeta, nBeta);

end
%% 6. REPORT RISULTATI

fprintf('\n============================================\n');
fprintf('       RISULTATI MONTE CARLO LTI\n');
fprintf('       %d controllori x %d campioni\n', Nc, N_campioni);
fprintf('============================================\n');

fprintf('Configurazioni totali: %d\n', Ntot);
fprintf('Configurazioni stabili: %d/%d\n', sum(Stable), Ntot);

% Indici utili
pitchSettled = Stable & PitchSettledWithinHorizon;
yawSettled   = Stable & YawSettledWithinHorizon;

%% ---------------------------------------------------------------
% RISULTATI PER CONTROLLATORE
% ---------------------------------------------------------------

for ic = 1:Nc

    idx = strcmp(resultController, string(controllerNames{ic}));

    fprintf('\n\n====================================================\n');
    fprintf('              CONTROLLATORE: %s\n', controllerNames{ic});
    fprintf('====================================================\n');

    %% STABILITÀ
    nStable = sum(Stable(idx));
    nTotCtrl = sum(idx);

    fprintf('\n--- STABILITA'' ---\n');
    fprintf('Configurazioni stabili = %d/%d\n', nStable, nTotCtrl);

    %% PITCH
    fprintf('\n--- PRESTAZIONI PITCH (ALPHA) ---\n');

    % Bandwidth
    validBW_pitch = idx & Stable & isfinite(PitchBW);

    if any(validBW_pitch)
        fprintf('Banda Passante media      = %.3f rad/s\n', ...
            mean(PitchBW(validBW_pitch)));
        fprintf('Banda Passante minima     = %.3f rad/s\n', ...
            min(PitchBW(validBW_pitch)));
        fprintf('Banda Passante massima    = %.3f rad/s\n', ...
            max(PitchBW(validBW_pitch)));
    else
        fprintf('Banda Passante media      = NaN rad/s\n');
    end

    % Settling
    validSettling_pitch = idx & pitchSettled;

    if any(validSettling_pitch)
        fprintf('Tempo Assestamento medio  = %.3f s\n', ...
            mean(PitchSettling(validSettling_pitch)));
        fprintf('Tempo Assestamento max     = %.3f s\n', ...
            max(PitchSettling(validSettling_pitch)));
        fprintf('Tempo Assestamento min     = %.3f s\n', ...
            min(PitchSettling(validSettling_pitch)));
    else
        fprintf('Tempo Assestamento max     = NaN s\n');
    end

    % Errore regime
    validSS_pitch = idx & Stable & isfinite(PitchSSerror);

    if any(validSS_pitch)
        fprintf('Errore a regime medio      = %.4f %%\n', ...
            mean(PitchSSerror(validSS_pitch)));
        fprintf('Errore a regime massimo    = %.4f %%\n', ...
            max(PitchSSerror(validSS_pitch)));
    else
        fprintf('Errore a regime massimo    = NaN %%\n');
    end

    % Disturbo
    validDist_pitch = idx & Stable & isfinite(PitchDistPeak);

    if any(validDist_pitch)
        fprintf('Picco massimo da Disturbo  = %.4e rad\n', ...
            max(PitchDistPeak(validDist_pitch)));
        fprintf('Picco medio da Disturbo    = %.4e rad\n', ...
            mean(PitchDistPeak(validDist_pitch)));
    else
        fprintf('Picco massimo da Disturbo  = NaN rad\n');
    end

    %% YAW
    fprintf('\n--- PRESTAZIONI YAW (BETA) ---\n');

    % Bandwidth
    validBW_yaw = idx & Stable & isfinite(YawBW);

    if any(validBW_yaw)
        fprintf('Banda Passante media      = %.3f rad/s\n', ...
            mean(YawBW(validBW_yaw)));
        fprintf('Banda Passante minima     = %.3f rad/s\n', ...
            min(YawBW(validBW_yaw)));
        fprintf('Banda Passante massima    = %.3f rad/s\n', ...
            max(YawBW(validBW_yaw)));
    else
        fprintf('Banda Passante media      = NaN rad/s\n');
    end

    % Settling
    validSettling_yaw = idx & yawSettled;

    if any(validSettling_yaw)
        fprintf('Tempo Assestamento medio  = %.3f s\n', ...
            mean(YawSettling(validSettling_yaw)));
        fprintf('Tempo Assestamento max     = %.3f s\n', ...
            max(YawSettling(validSettling_yaw)));
        fprintf('Tempo Assestamento min     = %.3f s\n', ...
            min(YawSettling(validSettling_yaw)));
    else
        fprintf('Tempo Assestamento max     = NaN s\n');
    end

    % Errore regime
    validSS_yaw = idx & Stable & isfinite(YawSSerror);

    if any(validSS_yaw)
        fprintf('Errore a regime medio      = %.4f %%\n', ...
            mean(YawSSerror(validSS_yaw)));
        fprintf('Errore a regime massimo    = %.4f %%\n', ...
            max(YawSSerror(validSS_yaw)));
    else
        fprintf('Errore a regime massimo    = NaN %%\n');
    end

    % Disturbo
    validDist_yaw = idx & Stable & isfinite(YawDistPeak);

    if any(validDist_yaw)
        fprintf('Picco massimo da Disturbo  = %.4e rad\n', ...
            max(YawDistPeak(validDist_yaw)));
        fprintf('Picco medio da Disturbo    = %.4e rad\n', ...
            mean(YawDistPeak(validDist_yaw)));
    else
        fprintf('Picco massimo da Disturbo  = NaN rad\n');
    end

    %% RECUPERO DAI DISTURBI

    fprintf('\n--- RECUPERO DAI DISTURBI ---\n');

    % Numero di campioni che rientrano stabilmente nella banda
    nRecAlpha = sum(PitchRecovered1s(idx & Stable));
    nRecBeta  = sum(YawRecovered1s(idx & Stable));

    fprintf('Pitch: rientro stabile nella banda per >= 1 s: %d/%d\n', ...
        nRecAlpha, nStable);

    fprintf('Yaw: rientro stabile nella banda per >= 1 s: %d/%d\n', ...
        nRecBeta, nStable);

    % Rispetto dei limiti temporali
    fprintf('Pitch: recupero entro 2 s: %d/%d\n', ...
        sum(PitchRecoveryPass(idx & Stable)), nStable);

    fprintf('Yaw: recupero entro 3 s: %d/%d\n', ...
        sum(YawRecoveryPass(idx & Stable)), nStable);

    % Statistiche dei tempi di recupero dei campioni recuperati
    validRecAlpha = idx & Stable & PitchRecovered1s;
    validRecBeta  = idx & Stable & YawRecovered1s;

    if any(validRecAlpha)
        fprintf('Trec pitch medio = %.3f s\n', ...
            mean(PitchRecoveryTime(validRecAlpha)));
        fprintf('Trec pitch massimo = %.3f s\n', ...
            max(PitchRecoveryTime(validRecAlpha)));
    else
        fprintf('Trec pitch: nessun recupero trovato entro 30 s.\n');
    end

    if any(validRecBeta)
        fprintf('Trec yaw medio = %.3f s\n', ...
            mean(YawRecoveryTime(validRecBeta)));
        fprintf('Trec yaw massimo = %.3f s\n', ...
            max(YawRecoveryTime(validRecBeta)));
    else
        fprintf('Trec yaw: nessun recupero trovato entro 30 s.\n');
    end

    % Errore residuo negli ultimi 5 secondi
    validResidAlpha = idx & Stable & isfinite(PitchResidualPeakLast5s);
    validResidBeta  = idx & Stable & isfinite(YawResidualPeakLast5s);

    if any(validResidAlpha)
        fprintf('Errore pitch massimo ultimi 5 s = %.4f deg\n', ...
            rad2deg(max(PitchResidualPeakLast5s(validResidAlpha))));
        fprintf('Errore pitch RMS medio ultimi 5 s = %.4f deg\n', ...
            rad2deg(mean(PitchResidualRMSLast5s(validResidAlpha))));
    end

    if any(validResidBeta)
        fprintf('Errore yaw massimo ultimi 5 s = %.4f deg\n', ...
            rad2deg(max(YawResidualPeakLast5s(validResidBeta))));
        fprintf('Errore yaw RMS medio ultimi 5 s = %.4f deg\n', ...
            rad2deg(mean(YawResidualRMSLast5s(validResidBeta))));
    end

    % Conformità RMS ultimi 5 secondi
    validRMSAlpha = idx & Stable & ...
        isfinite(PitchResidualRMSLast5s);

    validRMSBeta = idx & Stable & ...
        isfinite(YawResidualRMSLast5s);

    nRMSAlpha = sum(validRMSAlpha);
    nRMSBeta  = sum(validRMSBeta);

    nPassRMSAlpha = sum(PitchResidualRMSPass(validRMSAlpha));
    nPassRMSBeta  = sum(YawResidualRMSPass(validRMSBeta));

    fprintf('Pitch: RMS <= 0.3 deg: %d/%d (%.1f%%)\n', ...
        nPassRMSAlpha, nRMSAlpha, ...
        100*nPassRMSAlpha/max(nRMSAlpha,1));

    fprintf('Yaw: RMS <= 0.7 deg: %d/%d (%.1f%%)\n', ...
        nPassRMSBeta, nRMSBeta, ...
        100*nPassRMSBeta/max(nRMSBeta,1));
    if any(validRMSAlpha)
        fprintf('RMS pitch massimo tra i campioni = %.4f deg\n', ...
            rad2deg(max(PitchResidualRMSLast5s(validRMSAlpha))));
    end

    if any(validRMSBeta)
        fprintf('RMS yaw massimo tra i campioni = %.4f deg\n', ...
            rad2deg(max(YawResidualRMSLast5s(validRMSBeta))));
    end

    %% ACCOPPIAMENTO
    fprintf('\n--- ACCOPPIAMENTO ---\n');

    validCross1 = idx & Stable & isfinite(YawFromPitchPeak);
    validCross2 = idx & Stable & isfinite(PitchFromYawPeak);

    if any(validCross1)
        fprintf('Sbandamento Beta da Pitch  = %.4e rad\n', ...
            max(YawFromPitchPeak(validCross1)));
    else
        fprintf('Sbandamento Beta da Pitch  = NaN rad\n');
    end

    if any(validCross2)
        fprintf('Sbandamento Alpha da Yaw   = %.4e rad\n', ...
            max(PitchFromYawPeak(validCross2)));
    else
        fprintf('Sbandamento Alpha da Yaw   = NaN rad\n');
    end

    %% SFORZO DI CONTROLLO
    fprintf('\n--- SFORZO DI CONTROLLO ---\n');

    validU1 = idx & Stable & isfinite(max_u1_cmd);
    validU2 = idx & Stable & isfinite(max_u2_cmd);

    if any(validU1)
        worst_u1 = max(max_u1_cmd(validU1));
    else
        worst_u1 = NaN;
    end

    if any(validU2)
        worst_u2 = max(max_u2_cmd(validU2));
    else
        worst_u2 = NaN;
    end

    fprintf('Worst max |u1_cmd| = %.4f N (Limite = %.2f)\n', ...
        worst_u1, umax);

    fprintf('Worst max |u2_cmd| = %.4f N (Limite = %.2f)\n', ...
        worst_u2, umax);

    if any(idx & Stable)

        valid = idx & Stable;

        n_sat_ctrl = sum( ...
            min_u1_cmd(valid) < u1_min_lim | ...
            max_u1_signed(valid) > u1_max_lim | ...
            min_u2_cmd(valid) < u2_min_lim | ...
            max_u2_signed(valid) > u2_max_lim);

    else
        n_sat_ctrl = 0;
    end

    fprintf('Campioni oltre i limiti attuatori: %d/%d\n', ...
        n_sat_ctrl,nStable);

end


%% ---------------------------------------------------------------
% RIEPILOGO COMPLESSIVO
% ---------------------------------------------------------------

fprintf('\n\n============================================\n');
fprintf('         RIEPILOGO COMPLESSIVO\n');
fprintf('============================================\n');

fprintf('Configurazioni stabili totali: %d/%d\n', ...
    sum(Stable), Ntot);

validBW_pitch_all = Stable & isfinite(PitchBW);
validBW_yaw_all   = Stable & isfinite(YawBW);

if any(validBW_pitch_all)
    fprintf('Bandwidth pitch media      = %.3f rad/s\n', ...
        mean(PitchBW(validBW_pitch_all)));
end

if any(validBW_yaw_all)
    fprintf('Bandwidth yaw media        = %.3f rad/s\n', ...
        mean(YawBW(validBW_yaw_all)));
end

if any(pitchSettled)
    fprintf('Max settling pitch         = %.3f s\n', ...
        max(PitchSettling(pitchSettled)));
end

if any(yawSettled)
    fprintf('Max settling yaw           = %.3f s\n', ...
        max(YawSettling(yawSettled)));
end

fprintf('Max errore regime pitch    = %.4f %%\n', ...
    max(PitchSSerror(Stable)));

fprintf('Max errore regime yaw      = %.4f %%\n', ...
    max(YawSSerror(Stable)));

validOS_pitch = Stable & isfinite(PitchOvershoot);

if any(validOS_pitch)
    fprintf('Overshoot pitch medio           = %.4f %%\n', ...
        mean(PitchOvershoot(validOS_pitch)));
    fprintf('Overshoot pitch massimo         = %.4f %%\n', ...
        max(PitchOvershoot(validOS_pitch)));
else
    fprintf('Overshoot pitch massimo         = NaN %%\n');
end

validOS_yaw = Stable & isfinite(YawOvershoot);

if any(validOS_yaw)
    fprintf('Overshoot yaw medio           = %.4f %%\n', ...
        mean(YawOvershoot(validOS_yaw)));
    fprintf('Overshoot yaw massimo         = %.4f %%\n', ...
        max(YawOvershoot(validOS_yaw)));
else
    fprintf('Overshoot yaw massimo         = NaN %%\n');
end

fprintf('Max disturbo pitch         = %.4e rad\n', ...
    max(PitchDistPeak(Stable)));

fprintf('Max disturbo yaw           = %.4e rad\n', ...
    max(YawDistPeak(Stable)));

fprintf('Max beta da pitch          = %.4e rad\n', ...
    max(YawFromPitchPeak(Stable)));

fprintf('Max alpha da yaw            = %.4e rad\n', ...
    max(PitchFromYawPeak(Stable)));

fprintf('Worst max |u1_cmd|          = %.4f N\n', ...
    max(max_u1_cmd(Stable)));

fprintf('Worst max |u2_cmd|          = %.4f N\n', ...
    max(max_u2_cmd(Stable)));

%% Conteggio complessivo delle saturazioni
valid_sat = Stable;

n_sat = sum( ...
    min_u1_cmd(valid_sat) < u1_min_lim | ...
    max_u1_signed(valid_sat) > u1_max_lim | ...
    min_u2_cmd(valid_sat) < u2_min_lim | ...
    max_u2_signed(valid_sat) > u2_max_lim);

fprintf('Campioni totali in saturazione: %d/%d\n', ...
    n_sat, sum(valid_sat));
fprintf('Campioni totali in saturazione: %d/%d\n', ...
    n_sat, sum(Stable));

%% 7. PLOT CONFRONTO TRA CONTROLLORI
% ==========================================

ControllerPlot = categorical(resultController);

% Ordine dei controllori
labels = {'mixsyn','hinfsyn','PID+comp','mu-synthesis'};

%% Calcolo statistiche per controllore

meanPitchSettling = nan(Nc,1);
maxPitchSettling  = nan(Nc,1);

meanYawSettling = nan(Nc,1);
maxYawSettling  = nan(Nc,1);

meanPitchBW = nan(Nc,1);
meanYawBW   = nan(Nc,1);

meanPitchOS = nan(Nc,1);
maxPitchOS  = nan(Nc,1);

meanYawOS = nan(Nc,1);
maxYawOS  = nan(Nc,1);

meanYawFromPitch = nan(Nc,1);
maxYawFromPitch  = nan(Nc,1);

meanPitchFromYaw = nan(Nc,1);
maxPitchFromYaw  = nan(Nc,1);

meanPitchSSE = nan(Nc,1);
maxPitchSSE  = nan(Nc,1);

meanYawSSE = nan(Nc,1);
maxYawSSE  = nan(Nc,1);


for ic = 1:Nc

    idx = strcmp(resultController, string(controllerNames{ic})) & Stable;

    % --- Settling alpha ---
    v = PitchSettling(idx & isfinite(PitchSettling));
    if ~isempty(v)
        meanPitchSettling(ic) = mean(v);
        maxPitchSettling(ic)  = max(v);
    end

    % --- Settling beta ---
    v = YawSettling(idx & isfinite(YawSettling));
    if ~isempty(v)
        meanYawSettling(ic) = mean(v);
        maxYawSettling(ic)  = max(v);
    end

    % --- Bandwidth alpha ---
    v = PitchBW(idx & isfinite(PitchBW));
    if ~isempty(v)
        meanPitchBW(ic) = mean(v);
    end

    % --- Bandwidth beta ---
    v = YawBW(idx & isfinite(YawBW));
    if ~isempty(v)
        meanYawBW(ic) = mean(v);
    end

    % --- Overshoot alpha ---
    v = PitchOvershoot(idx & isfinite(PitchOvershoot));
    if ~isempty(v)
        meanPitchOS(ic) = mean(v);
        maxPitchOS(ic)  = max(v);
    end

    % --- Overshoot beta ---
    v = YawOvershoot(idx & isfinite(YawOvershoot));
    if ~isempty(v)
        meanYawOS(ic) = mean(v);
        maxYawOS(ic)  = max(v);
    end

    % --- Cross coupling alpha -> beta ---
    v = YawFromPitchPeak(idx & isfinite(YawFromPitchPeak));
    if ~isempty(v)
        meanYawFromPitch(ic) = mean(v);
        maxYawFromPitch(ic)  = max(v);
    end

    % --- Cross coupling beta -> alpha ---
    v = PitchFromYawPeak(idx & isfinite(PitchFromYawPeak));
    if ~isempty(v)
        meanPitchFromYaw(ic) = mean(v);
        maxPitchFromYaw(ic)  = max(v);
    end

    % --- Errore regime alpha ---
    v = PitchSSerror(idx & isfinite(PitchSSerror));
    if ~isempty(v)
        meanPitchSSE(ic) = mean(v);
        maxPitchSSE(ic)  = max(v);
    end

    % --- Errore regime beta ---
    v = YawSSerror(idx & isfinite(YawSSerror));
    if ~isempty(v)
        meanYawSSE(ic) = mean(v);
        maxYawSSE(ic)  = max(v);
    end

end


%% ==========================================
% FIGURA 1 - SETTLING TIME
% ==========================================

figure('Name','Settling Time','Color','w');

subplot(1,2,1)

bar(meanPitchSettling)
hold on
plot(1:Nc,maxPitchSettling,'k^','MarkerSize',7,'LineWidth',1.2)
hold off

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Settling time \alpha [s]')
title('Settling time \alpha')
legend('Media','Massimo','Location','best')
grid on


subplot(1,2,2)

bar(meanYawSettling)
hold on
plot(1:Nc,maxYawSettling,'k^','MarkerSize',7,'LineWidth',1.2)
hold off

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Settling time \beta [s]')
title('Settling time \beta')
legend('Media','Massimo','Location','best')
grid on


%% ==========================================
% FIGURA 2 - BANDWIDTH
% ==========================================

figure('Name','Bandwidth','Color','w');

subplot(1,2,1)

bar(meanPitchBW)

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Bandwidth \alpha [rad/s]')
title('Bandwidth \alpha')
grid on


subplot(1,2,2)

bar(meanYawBW)

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Bandwidth \beta [rad/s]')
title('Bandwidth \beta')
grid on


%% ==========================================
% FIGURA 3 - OVERSHOOT
% ==========================================

figure('Name','Overshoot','Color','w');

subplot(1,2,1)

bar(meanPitchOS)
hold on
plot(1:Nc,maxPitchOS,'k^','MarkerSize',7,'LineWidth',1.2)
hold off

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Overshoot \alpha [%]')
title('Overshoot \alpha')
legend('Media','Massimo','Location','best')
grid on


subplot(1,2,2)

bar(meanYawOS)
hold on
plot(1:Nc,maxYawOS,'k^','MarkerSize',7,'LineWidth',1.2)
hold off

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Overshoot \beta [%]')
title('Overshoot \beta')
legend('Media','Massimo','Location','best')
grid on


%% ==========================================
% FIGURA 4 - CROSS COUPLING
% ==========================================

figure('Name','Cross Coupling','Color','w');

subplot(1,2,1)

bar(meanYawFromPitch)
hold on
plot(1:Nc,maxYawFromPitch,'k^','MarkerSize',7,'LineWidth',1.2)
hold off

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('\beta \leftarrow \alpha [rad]')
title('\alpha \rightarrow \beta')
legend('Media','Massimo','Location','best')
grid on


subplot(1,2,2)

bar(meanPitchFromYaw)
hold on
plot(1:Nc,maxPitchFromYaw,'k^','MarkerSize',7,'LineWidth',1.2)
hold off

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('\alpha \leftarrow \beta [rad]')
title('\beta \rightarrow \alpha')
legend('Media','Massimo','Location','best')
grid on

%% ==========================================
% FIGURA 5 - ERRORE A REGIME
% ==========================================

figure('Name','Errore a regime','Color','w');

% --- Alpha ---
subplot(1,2,1)

bar(meanPitchSSE)

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Errore a regime \alpha [%]')
title('Errore a regime \alpha')
grid on

for ic = 1:Nc
    if isfinite(meanPitchSSE(ic))
        text(ic, meanPitchSSE(ic), ...
            sprintf('%.3g %%',meanPitchSSE(ic)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom');
    end
end


% --- Beta ---
subplot(1,2,2)

bar(meanYawSSE)

set(gca,'XTick',1:Nc,'XTickLabel',labels)
ylabel('Errore a regime \beta [%]')
title('Errore a regime \beta')
grid on

for ic = 1:Nc
    if isfinite(meanYawSSE(ic))
        text(ic, meanYawSSE(ic), ...
            sprintf('%.3g %%',meanYawSSE(ic)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom');
    end
end

sgtitle('Errore a regime - Monte Carlo');
MonteCarloResults = table( ...
    resultController, ...
    resultSample, ...
    Stable, ...
    PitchSettling, ...
    YawSettling, ...
    PitchSettledWithinHorizon, ...
    YawSettledWithinHorizon, ...
    PitchBW, ...
    YawBW, ...
    PitchOvershoot, ...
    YawOvershoot, ...
    PitchSSerror, ...
    YawSSerror, ...
    YawFromPitchPeak, ...
    PitchFromYawPeak, ...
    PitchDistPeak, ...
    YawDistPeak, ...
    PitchRecoveryTime, ...
    YawRecoveryTime, ...
    PitchRecovered1s, ...
    YawRecovered1s, ...
    PitchRecoveryPass, ...
    YawRecoveryPass, ...
    PitchResidualPeakLast5s, ...
    YawResidualPeakLast5s, ...
    PitchResidualRMSLast5s, ...
    YawResidualRMSLast5s, ...
    PitchResidualRMSPass, ...
    YawResidualRMSPass, ...
    'VariableNames',{ ...
    'Controller'
    'Sample'
    'Stable'
    'PitchSettling'
    'YawSettling'
    'PitchSettledWithinHorizon'
    'YawSettledWithinHorizon'
    'PitchBandwidth'
    'YawBandwidth'
    'PitchOvershootPercent'
    'YawOvershootPercent'
    'PitchSteadyStateErrorPercent'
    'YawSteadyStateErrorPercent'
    'YawFromPitchPeak'
    'PitchFromYawPeak'
    'PitchDisturbancePeak'
    'YawDisturbancePeak'
    'PitchRecoveryTime_s'
    'YawRecoveryTime_s'
    'PitchRecoveredWithinBand1s'
    'YawRecoveredWithinBand1s'
    'PitchRecoveryWithinSpec'
    'YawRecoveryWithinSpec'
    'PitchResidualPeakLast5s_rad'
    'YawResidualPeakLast5s_rad'
    'PitchResidualRMSLast5s_rad'
    'YawResidualRMSLast5s_rad'
    'PitchResidualRMSPass'
    'YawResidualRMSPass'
    });

disp(MonteCarloResults);

writetable( ...
    MonteCarloResults, ...
    'MonteCarloRobustness.csv');

function [ts,settled] = settlingTimeFromTrace(t,y,yss,threshold)

    t = t(:);
    y = y(:);

    scale = max(abs(yss),1e-9);
    band = threshold*scale;
    outside = abs(y-yss) > band;

    lastOutside = find(outside,1,'last');

    if isempty(lastOutside)
        ts = t(1);
        settled = true;
    elseif lastOutside < numel(t)
        ts = t(lastOutside+1);
        settled = true;
    else
        ts = NaN;
        settled = false;
    end
end

function [trec, recovered] = recoveryTimeFromError( ...
    t, e, td, band, dwell)

t = t(:);
e = e(:);

trec = NaN;
recovered = false;

if numel(t) < 2
    return;
end

dt = median(diff(t));
nHold = ceil(dwell/dt) + 1;

k0 = find(t >= td, 1, 'first');

if isempty(k0)
    return;
end

insideBand = abs(e) <= band;

% Cerca la prima uscita dalla banda dopo l'applicazione
kExitRel = find(~insideBand(k0:end), 1, 'first');

% Se l'errore non esce mai dalla banda, non c'è
% una fase di recupero da misurare: la banda è rispettata.
if isempty(kExitRel)
    trec = 0;
    recovered = true;
    return;
end

kExit = k0 + kExitRel - 1;
kStartSearch = kExit + 1;

if kStartSearch + nHold - 1 > numel(t)
    return;
end

% Cerca il primo rientro seguito da almeno dwell secondi
for k = kStartSearch:(numel(t) - nHold + 1)

    idxHold = k:(k + nHold - 1);

    if all(insideBand(idxHold))
        trec = t(k) - td;
        recovered = true;
        return;
    end
end
end