%% ========================================================================
% MU-SYNTHESIS DIRETTA NON-LUMPED (COLD-START) - ELICOTTERO 2-DOF
%
% Questo script esegue SOLTANTO la sintesi diretta sul modello non-lumped.
%
% Il controllore viene salvato in un file MAT separato, senza sovrascrivere
% altri risultati precedenti. Se il file di output dedicato esiste gia',
% viene creato un nuovo file con timestamp.
% ========================================================================

close all;
clc;

fprintf('============================================================\n');
fprintf(' MU-SYNTHESIS DIRETTA NON-LUMPED - COLD START\n');
fprintf('============================================================\n');

%% 1. CARICAMENTO DEL SETUP COMUNE (NESSUN FILE LUMPED)

if ~isfile('HINF_workspace.mat')
    error(['File HINF_workspace.mat non trovato nella cartella corrente. ' ...
           'Portare lo script nella cartella del progetto o impostare il percorso corretto.']);
end

load('HINF_workspace.mat');

requiredVars = {'G_unc_mu_scaled','G_unc_mu','WS','WU','WT', ...
                'Du','Dy_inv'};
for kVar = 1:numel(requiredVars)
    if ~exist(requiredVars{kVar},'var')
        error('La variabile %s manca in HINF_workspace.mat. Rieseguire HINF_SETUP.', ...
              requiredVars{kVar});
    end
end

requiredUncertainty = {'J_alpha','l','omega_n1','omega_n2', ...
                       'tau_d1','tau_d2','Delta_Pade1','Delta_Pade2'};
availableUncertainty = fieldnames(G_unc_mu_scaled.Uncertainty);
missingUncertainty = setdiff(requiredUncertainty,availableUncertainty);
if ~isempty(missingUncertainty)
    error(['G_unc_mu_scaled non contiene tutte le incertezze attese. ' ...
           'Rieseguire HINF_SETUP. Mancano: %s'], ...
           strjoin(missingUncertainty,', '));
end

fprintf('Caricato HINF_workspace.mat.\n');
fprintf('La sintesi lumped e'' disattivata: non verra'' chiamata musyn sul plant lumped.\n');

%% 2. IMPIANTO GENERALIZZATO NON-LUMPED

Gmu_nonlumped_synth = G_unc_mu_scaled;
Gmu_nonlumped_synth.InputName  = {'ubar1','ubar2'};
Gmu_nonlumped_synth.OutputName = {'alpha_bar','beta_bar'};

% Nomi dei segnali dei pesi, coerenti con lo script originale.
WS = ss(WS);
WS.InputName  = {'e1','e2'};
WS.OutputName = {'zS1','zS2'};

WU = ss(WU);
WU.InputName  = {'ubar1','ubar2'};
WU.OutputName = {'zU1','zU2'};

WT = ss(WT);
WT.InputName  = {'alpha_bar','beta_bar'};
WT.OutputName = {'zT1','zT2'};

SumE1 = sumblk('e1 = r1 - alpha_bar');
SumE2 = sumblk('e2 = r2 - beta_bar');

P_mu_nonlumped = connect( ...
    Gmu_nonlumped_synth, ...
    WS, WU, WT, ...
    SumE1, SumE2, ...
    {'r1','r2','ubar1','ubar2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2','e1','e2'});

fprintf('Plant generalizzato non-lumped creato.\n');

%% 3. OPZIONI MUSYN (IDENTICHE AL RAMO NON-LUMPED DELLO SCRIPT ORIGINALE)

nmeas = 2;
ncont = 2;

optsMU_nonlumped = musynOptions( ...
    'Display','short', ...
    'MixedMU','on', ...
    'FullDG',false, ...
    'FitOrder',[3 2], ...
    'MaxIter',40, ...
    'TolPerf',0.01);

% Output dedicato: non tocca mai i MAT lumped o il file MU_controller.mat.
outputFile = 'MU_controller_NONLUMPED_COLDSTART_ONLY.mat';
if isfile(outputFile)
    [outDir,outName,outExt] = fileparts(outputFile);
    if isempty(outDir)
        outDir = '.';
    end
    outputFile = fullfile(outDir, ...
        [outName '_' datestr(now,'yyyymmdd_HHMMSS') outExt]);
end

%% 4. SINTESI DIRETTA NON-LUMPED DA COLD-START

fprintf('\n============================================================\n');
fprintf('AVVIO MUSYN SUL SOLO MODELLO NON-LUMPED\n');
fprintf('Cold-start: nessun controllore iniziale passato a musyn.\n');
fprintf('MaxIter = 8 | FullDG = false | FitOrder = [3 2] | TolPerf = 0.01\n');
fprintf('============================================================\n');

rng(10);
[K_mu_scaled_nonlumped, ...
 CLperf_mu_nonlumped, ...
 info_mu_nonlumped] = musyn( ...
    P_mu_nonlumped,nmeas,ncont,optsMU_nonlumped);

K_mu_scaled_nonlumped = ss(K_mu_scaled_nonlumped);

fprintf('\nSintesi non-lumped conclusa.\n');
fprintf('Ordine K non-lumped = %d\n',order(K_mu_scaled_nonlumped));
fprintf('Performance restituita da musyn = %.8g\n',CLperf_mu_nonlumped);
fprintf('Controllore stabile = %d\n',isstable(K_mu_scaled_nonlumped));

%% 5. DESCALATURA NELLE COORDINATE FISICHE

% Coerente con lo script originale: K_phys = Du * K_scaled * Dy_inv.
K_mu_nonlumped_phys = Du * K_mu_scaled_nonlumped * Dy_inv;
K_mu_nonlumped_phys.InputName  = {'e1','e2'};
K_mu_nonlumped_phys.OutputName = {'u1','u2'};

%% 6. SALVATAGGIO IMMEDIATO DEL CONTROLLore
% Il salvataggio avviene prima delle verifiche: il controllore resta salvato
% anche se una successiva analisi di robustezza dovesse generare un errore.

save(outputFile, ...
    'K_mu_scaled_nonlumped', ...
    'K_mu_nonlumped_phys', ...
    'CLperf_mu_nonlumped', ...
    'info_mu_nonlumped', ...
    'optsMU_nonlumped', ...
    'P_mu_nonlumped', ...
    'nmeas','ncont','Du','Dy_inv');

fprintf('\nControllore salvato subito in:\n%s\n',outputFile);
fprintf('Nessun file lumped o MAT precedente e'' stato sovrascritto.\n');

%% 7. VERIFICHE SUL MODELLO NON-LUMPED
% Le verifiche sono separate dalla sintesi. In caso di errore, il controllore
% e' gia' stato salvato nel file dedicato.

optsWC = wcOptions('Display','off','MussvOptions','a');

% Closed-loop pesato sul modello non-lumped normalizzato.
Gmu_nonlumped = G_unc_mu_scaled;
Gmu_nonlumped.InputName  = {'u1','u2'};
Gmu_nonlumped.OutputName = {'y1','y2'};

Ktest = ss(K_mu_scaled_nonlumped);
Ktest.InputName  = {'e1','e2'};
Ktest.OutputName = {'u1','u2'};

W_S = ss(WS);
W_S.InputName  = {'e1','e2'};
W_S.OutputName = {'zS1','zS2'};

W_U = ss(WU);
W_U.InputName  = {'u1','u2'};
W_U.OutputName = {'zU1','zU2'};

W_T = ss(WT);
W_T.InputName  = {'y1','y2'};
W_T.OutputName = {'zT1','zT2'};

sum1 = sumblk('e1 = r1 - y1');
sum2 = sumblk('e2 = r2 - y2');

CL_nonlumped = connect( ...
    Gmu_nonlumped,Ktest,W_S,W_U,W_T,sum1,sum2, ...
    {'r1','r2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2'});

% NS e NP nominali.
CLnom_nonlumped = getNominal(CL_nonlumped);
NS_nonlumped = isstable(CLnom_nonlumped);
if NS_nonlumped
    gammaNom_nonlumped = hinfnorm(minreal(CLnom_nonlumped,1e-5));
else
    gammaNom_nonlumped = Inf;
end
NP_nonlumped = gammaNom_nonlumped < 1;

fprintf('\n============================================================\n');
fprintf('RISULTATI NOMINALI SUL NON-LUMPED\n');
fprintf('============================================================\n');
fprintf('NS nominale = %d\n',NS_nonlumped);
fprintf('Gamma NP nominale = %.8g\n',gammaNom_nonlumped);
fprintf('NP < 1 = %d\n',NP_nonlumped);

% RP robusta. Il controllore e' gia' stato salvato se wcgain fallisce.
try
    [wc_nonlumped_synth,wcu_nonlumped_synth] = wcgain(CL_nonlumped,optsWC);
    RP_nonlumped_lower = wc_nonlumped_synth.LowerBound;
    RP_nonlumped_upper = wc_nonlumped_synth.UpperBound;

    fprintf('\n============================================================\n');
    fprintf('ROBUST PERFORMANCE SUL NON-LUMPED\n');
    fprintf('============================================================\n');
    fprintf('RP lower = %.8g\n',RP_nonlumped_lower);
    fprintf('RP upper = %.8g\n',RP_nonlumped_upper);
    fprintf('RP certificata < 1 = %d\n',RP_nonlumped_upper < 1);

    save(outputFile, ...
        'CL_nonlumped','CLnom_nonlumped', ...
        'NS_nonlumped','gammaNom_nonlumped','NP_nonlumped', ...
        'wc_nonlumped_synth','wcu_nonlumped_synth', ...
        'RP_nonlumped_lower','RP_nonlumped_upper','-append');
catch ME
    warning('Verifica wcgain non completata. Il controllore e'' comunque salvato. Dettaglio: %s',ME.message);
end

% RS sul modello fisico non-lumped con controllore descalato.
try
    Gmu_phys = G_unc_mu;
    Gmu_phys.InputName  = {'u1','u2'};
    Gmu_phys.OutputName = {'alpha','beta'};

    SumE1_phys = sumblk('e1 = r1 - alpha');
    SumE2_phys = sumblk('e2 = r2 - beta');

    CL_unc_mu = connect( ...
        Gmu_phys,K_mu_nonlumped_phys,SumE1_phys,SumE2_phys, ...
        {'r1','r2'},{'alpha','beta'});

    [stabmarg_nonlumped,~,info_RS_nonlumped] = robuststab(CL_unc_mu);

    fprintf('\n============================================================\n');
    fprintf('ROBUST STABILITY SUL NON-LUMPED\n');
    fprintf('============================================================\n');
    fprintf('RS lower = %.8g\n',stabmarg_nonlumped.LowerBound);
    fprintf('RS upper = %.8g\n',stabmarg_nonlumped.UpperBound);
    fprintf('RS certificata = %d\n',stabmarg_nonlumped.LowerBound > 1);

    save(outputFile, ...
        'CL_unc_mu','stabmarg_nonlumped','info_RS_nonlumped','-append');
catch ME
    warning('Verifica robuststab non completata. Il controllore e'' comunque salvato. Dettaglio: %s',ME.message);
end

fprintf('\n============================================================\n');
fprintf('FINE SCRIPT NON-LUMPED COLD-START\n');
fprintf('File risultato: %s\n',outputFile);
fprintf('============================================================\n');
