close all
clc

%% ==========================================
% ANALISI PRESTAZIONI LQG - SIMULINK
%
% tipo_LQG:
% 1 = LQG 1-DOF
% 2 = LQG 2-DOF senza integratore
% 3 = LQG 2-DOF con integratore
%
% tipo_test:
% 1 = tracking nominale, senza rumore/disturbo
% 2 = tracking + disturbo + rumore
% 3 = disturbo + rumore
%
% IMPORTANTE:
% La simulazione LINEARE e NON LINEARE devono essere
% già state eseguite e i relativi segnali devono essere
% presenti contemporaneamente in "out".
%% ==========================================

tipo_LQG  = 3;
tipo_test = 2;


%% ==========================================
% 1. SELEZIONE SEGNALI MODELLO LINEARE
%% ==========================================

switch tipo_LQG

    case 1

        alpha_lin_obj = out.alpha_LQG;
        beta_lin_obj  = out.beta_LQG;

        u_cmd_lin_obj = out.u_cmd_LQG;

        nome = 'LQG 1-DOF senza integratore';


    case 2

        alpha_lin_obj = out.alpha_LQG2;
        beta_lin_obj  = out.beta_LQG2;

        u_cmd_lin_obj = out.u_cmd_LQG;

        nome = 'LQG 2-DOF senza integratore';


    case 3

        alpha_lin_obj = out.alpha_LQG_int;
        beta_lin_obj  = out.beta_LQG_int;

        u_cmd_lin_obj = out.u_cmd_LQG;

        nome = 'LQG 2-DOF con integratore';


    otherwise

        error('tipo_LQG deve essere 1, 2 oppure 3.')

end


%% ==========================================
% 2. SELEZIONE SEGNALI MODELLO NON LINEARE
%% ==========================================

switch tipo_LQG

    case 1

        alpha_nl_obj = out.alpha_LQG_nl;
        beta_nl_obj  = out.beta_LQG_nl;

        u_cmd_nl_obj = out.u_cmd_LQG_nl;


    case 2

        alpha_nl_obj = out.alpha_LQG2_nl;
        beta_nl_obj  = out.beta_LQG2_nl;

        u_cmd_nl_obj = out.u_cmd_LQG2_nl;


    case 3

        alpha_nl_obj = out.alpha_LQG_int_nl;
        beta_nl_obj  = out.beta_LQG_int_nl;

        u_cmd_nl_obj = out.u_cmd_LQG_nl;


    otherwise

        error('tipo_LQG deve essere 1, 2 oppure 3.')

end


%% ==========================================
% 3. ESTRAZIONE VETTORI LINEARI
%% ==========================================

t_alpha_lin     = alpha_lin_obj.Time(:);
alpha_lin_value = alpha_lin_obj.Data(:);

t_beta_lin      = beta_lin_obj.Time(:);
beta_lin_value  = beta_lin_obj.Data(:);


%% ==========================================
% 4. ESTRAZIONE VETTORI NON LINEARI
%% ==========================================

t_alpha_nl     = alpha_nl_obj.Time(:);
alpha_nl_value = alpha_nl_obj.Data(:);

t_beta_nl      = beta_nl_obj.Time(:);
beta_nl_value  = beta_nl_obj.Data(:);


%% ==========================================
% 5. RIFERIMENTI
%% ==========================================

switch tipo_test

    case 1

        r_alpha = deg2rad(3);
        r_beta  = 0;


    case 2

        r_alpha = deg2rad(3);
        r_beta  = 0;


    case 3

        r_alpha = 0;
        r_beta  = 0;


    otherwise

        error('tipo_test deve essere 1, 2 oppure 3.')

end


%% ==========================================
% 6. ANALISI MODELLO LINEARE
%% ==========================================

fprintf('\n\n')
fprintf('##################################################\n')
fprintf('# ANALISI MODELLO LINEARE\n')
fprintf('##################################################\n')

ris_lin = calcola_prestazioni_LQG( ...
    alpha_lin_obj, ...
    beta_lin_obj, ...
    tipo_test, ...
    [nome ' - Lineare']);


%% ==========================================
% 7. ANALISI MODELLO NON LINEARE
%% ==========================================

fprintf('\n\n')
fprintf('##################################################\n')
fprintf('# ANALISI MODELLO NON LINEARE\n')
fprintf('##################################################\n')

ris_nl = calcola_prestazioni_LQG( ...
    alpha_nl_obj, ...
    beta_nl_obj, ...
    tipo_test, ...
    [nome ' - Non lineare']);


%% ==========================================
% 8. COMANDO DI CONTROLLO
%% ==========================================

t_u_lin = u_cmd_lin_obj.Time(:);

u_cmd_lin1 = u_cmd_lin_obj.Data(:,1);
u_cmd_lin2 = u_cmd_lin_obj.Data(:,2);


t_u_nl = u_cmd_nl_obj.Time(:);

u_cmd_nl1 = u_cmd_nl_obj.Data(:,1);
u_cmd_nl2 = u_cmd_nl_obj.Data(:,2);


%% ==========================================
% 9. DISTURBI AERODINAMICI
%% ==========================================

if isprop(out,'aereodynamic_disturbances') || ...
        isfield(out,'aereodynamic_disturbances')

    aero_obj = out.aereodynamic_disturbances;

    t_aero = aero_obj.Time(:);

    dist_aero_alpha = aero_obj.Data(:,1);
    dist_aero_beta  = aero_obj.Data(:,2);

else

    warning('Segnale aereodynamic_disturbances non trovato.')

    t_aero = [];
    dist_aero_alpha = [];
    dist_aero_beta = [];

end


%% ==========================================
% 10. ERRORE DI LINEARIZZAZIONE
%% ==========================================

% Portiamo il modello lineare sulla griglia temporale
% del modello non lineare.

alpha_lin_interp = interp1( ...
    t_alpha_lin, ...
    alpha_lin_value, ...
    t_alpha_nl, ...
    'linear', ...
    'extrap');


beta_lin_interp = interp1( ...
    t_beta_lin, ...
    beta_lin_value, ...
    t_beta_nl, ...
    'linear', ...
    'extrap');


err_lin_alpha = alpha_nl_value - alpha_lin_interp;

err_lin_beta = beta_nl_value - beta_lin_interp;


%% ==========================================
% 11. CONFRONTO DEI RISULTATI NUMERICI
%% ==========================================

fprintf('\n\n')
fprintf('==================================================\n')
fprintf(' CONFRONTO LINEARE vs NON LINEARE\n')
fprintf('==================================================\n')

fprintf('\n--- ERRORE RMS ---\n')

fprintf('RMS alpha LIN = %.6f rad\n', ...
    ris_lin.RMS_alpha);

fprintf('RMS alpha NL  = %.6f rad\n', ...
    ris_nl.RMS_alpha);

fprintf('RMS beta LIN  = %.6f rad\n', ...
    ris_lin.RMS_beta);

fprintf('RMS beta NL   = %.6f rad\n', ...
    ris_nl.RMS_beta);


fprintf('\n--- PICCO ---\n')

fprintf('Peak alpha LIN = %.6f rad\n', ...
    ris_lin.peak_alpha);

fprintf('Peak alpha NL  = %.6f rad\n', ...
    ris_nl.peak_alpha);

fprintf('Peak beta LIN  = %.6f rad\n', ...
    ris_lin.peak_beta);

fprintf('Peak beta NL   = %.6f rad\n', ...
    ris_nl.peak_beta);


fprintf('\n--- ERRORE FINALE ---\n')

fprintf('Final error alpha LIN = %.6f rad\n', ...
    ris_lin.final_error_alpha);

fprintf('Final error alpha NL  = %.6f rad\n', ...
    ris_nl.final_error_alpha);

fprintf('Final error beta LIN  = %.6f rad\n', ...
    ris_lin.final_error_beta);

fprintf('Final error beta NL   = %.6f rad\n', ...
    ris_nl.final_error_beta);


%% ==========================================
% 12. GRAFICO TRACKING LINEARE vs NON LINEARE
%% ==========================================

figure( ...
    'Name',[nome ' - Linear vs Nonlinear'], ...
    'Color','w')


subplot(2,1,1)

plot( ...
    t_alpha_lin, ...
    rad2deg(alpha_lin_value), ...
    '--', ...
    'LineWidth',1.5)

hold on

plot( ...
    t_alpha_nl, ...
    rad2deg(alpha_nl_value), ...
    '-', ...
    'LineWidth',1.5)

yline( ...
    rad2deg(r_alpha), ...
    ':', ...
    'LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('\alpha [deg]')

title('Pitch: linearized vs nonlinear model')

legend( ...
    'Linearized', ...
    'Nonlinear', ...
    'Reference', ...
    'Location','best')


subplot(2,1,2)

plot( ...
    t_beta_lin, ...
    rad2deg(beta_lin_value), ...
    '--', ...
    'LineWidth',1.5)

hold on

plot( ...
    t_beta_nl, ...
    rad2deg(beta_nl_value), ...
    '-', ...
    'LineWidth',1.5)

yline( ...
    rad2deg(r_beta), ...
    ':', ...
    'LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('\beta [deg]')

title('Yaw: linearized vs nonlinear model')

legend( ...
    'Linearized', ...
    'Nonlinear', ...
    'Reference', ...
    'Location','best')


%% ==========================================
% 13. ERRORE DI LINEARIZZAZIONE
%% ==========================================

figure( ...
    'Name',[nome ' - Linearization Error'], ...
    'Color','w')


subplot(2,1,1)

plot( ...
    t_alpha_nl, ...
    rad2deg(err_lin_alpha), ...
    'LineWidth',1.5)

hold on

yline(0,'k--','LineWidth',0.5)

grid on

xlabel('Time [s]')
ylabel('\alpha_{NL}-\alpha_{LIN} [deg]')

title('Pitch linearization error')


subplot(2,1,2)

plot( ...
    t_beta_nl, ...
    rad2deg(err_lin_beta), ...
    'LineWidth',1.5)

hold on

yline(0,'k--','LineWidth',0.5)

grid on

xlabel('Time [s]')
ylabel('\beta_{NL}-\beta_{LIN} [deg]')

title('Yaw linearization error')


%% ==========================================
% 14. SFORZO DI CONTROLLO
% ==========================================

figure( ...
    'Name',[nome ' - Control Effort Comparison'], ...
    'Color','w')


subplot(2,1,1)

plot( ...
    t_u_lin, ...
    u_cmd_lin1, ...
    '--', ...
    'LineWidth',1.5)

hold on

plot( ...
    t_u_nl, ...
    u_cmd_nl1, ...
    '-', ...
    'LineWidth',1.5)

grid on

title('Main-rotor control command: Linear vs Nonlinear')

ylabel('\Delta F_1 [N]')
xlabel('Time [s]')

legend( ...
    'Linearized', ...
    'Nonlinear', ...
    'Location','best')


subplot(2,1,2)

plot( ...
    t_u_lin, ...
    u_cmd_lin2, ...
    '--', ...
    'LineWidth',1.5)

hold on

plot( ...
    t_u_nl, ...
    u_cmd_nl2, ...
    '-', ...
    'LineWidth',1.5)

grid on

title('Tail-rotor control command: Linear vs Nonlinear')

ylabel('\Delta F_2 [N]')
xlabel('Time [s]')

legend( ...
    'Linearized', ...
    'Nonlinear', ...
    'Location','best')


%% ==========================================
% 15. COMANDO DI CONTROLLO LINEARE
%% ==========================================

figure( ...
    'Name',[nome ' - Control Effort Linear'], ...
    'Color','w')


subplot(2,1,1)

plot( ...
    t_u_lin, ...
    u_cmd_lin1, ...
    'LineWidth',1.5)

grid on

title('Main-rotor control command - Linear')

ylabel('\Delta F_1 [N]')
xlabel('Time [s]')


subplot(2,1,2)

plot( ...
    t_u_lin, ...
    u_cmd_lin2, ...
    'LineWidth',1.5)

grid on

title('Tail-rotor control command - Linear')

ylabel('\Delta F_2 [N]')
xlabel('Time [s]')


%% ==========================================
% 16. COMANDO DI CONTROLLO NON LINEARE
%% ==========================================

figure( ...
    'Name',[nome ' - Control Effort Nonlinear'], ...
    'Color','w')


subplot(2,1,1)

plot( ...
    t_u_nl, ...
    u_cmd_nl1, ...
    'LineWidth',1.5)

grid on

title('Main-rotor control command - Nonlinear')

ylabel('\Delta F_1 [N]')
xlabel('Time [s]')


subplot(2,1,2)

plot( ...
    t_u_nl, ...
    u_cmd_nl2, ...
    'LineWidth',1.5)

grid on

title('Tail-rotor control command - Nonlinear')

ylabel('\Delta F_2 [N]')
xlabel('Time [s]')


%% ==========================================
% 17. DISTURBI AERODINAMICI
%% ==========================================

if ~isempty(t_aero)

    figure( ...
        'Name',[nome ' - Aerodynamic Disturbances'], ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        t_aero, ...
        dist_aero_alpha, ...
        'LineWidth',1.5)

    grid on

    title('Pitch aerodynamic disturbance')

    ylabel('d_\alpha [N m]')
    xlabel('Time [s]')


    subplot(2,1,2)

    plot( ...
        t_aero, ...
        dist_aero_beta, ...
        'LineWidth',1.5)

    grid on

    title('Yaw aerodynamic disturbance')

    ylabel('d_\beta [N m]')
    xlabel('Time [s]')

end


%% ==========================================
% 18. TABELLA RIASSUNTIVA
%% ==========================================

Metriche = { ...
    'RMS errore alpha'; ...
    'RMS errore beta'; ...
    'Peak alpha'; ...
    'Peak beta'; ...
    'Errore finale alpha'; ...
    'Errore finale beta'; ...
    'Tempo recupero alpha'; ...
    'Tempo recupero beta'};

Lineare = [ ...
    ris_lin.RMS_alpha; ...
    ris_lin.RMS_beta; ...
    ris_lin.peak_alpha; ...
    ris_lin.peak_beta; ...
    ris_lin.final_error_alpha; ...
    ris_lin.final_error_beta; ...
    ris_lin.T_rec_alpha; ...
    ris_lin.T_rec_beta];

NonLineare = [ ...
    ris_nl.RMS_alpha; ...
    ris_nl.RMS_beta; ...
    ris_nl.peak_alpha; ...
    ris_nl.peak_beta; ...
    ris_nl.final_error_alpha; ...
    ris_nl.final_error_beta; ...
    ris_nl.T_rec_alpha; ...
    ris_nl.T_rec_beta];

Tabella_confronto = table( ...
    Metriche, ...
    Lineare, ...
    NonLineare);

disp(' ')
disp('==================================================')
disp(' TABELLA CONFRONTO')
disp('==================================================')
disp(Tabella_confronto)