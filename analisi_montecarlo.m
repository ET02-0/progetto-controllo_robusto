%% ==========================================
% MONTE CARLO LQG IN SIMULINK
% ==========================================

disp('--- Avvio Simulazione Monte Carlo per LQG in Simulink ---')

N_campioni = 50;
rng('default');

%% 1. SALVO IL MODELLO NOMINALE

A_nom_saved = A_nom;
B_nom_saved = B_nom;
Cy_saved   = Cy;
Dy_saved   = Dy;


% ==========================================
% METRICHE MONTE CARLO
% ==========================================

max_alpha     = zeros(N_campioni,1);
max_beta      = zeros(N_campioni,1);

rms_alpha     = zeros(N_campioni,1);
rms_beta      = zeros(N_campioni,1);

final_alpha   = zeros(N_campioni,1);
final_beta    = zeros(N_campioni,1);

beta_dev_ss = zeros(N_campioni,1);

T_rec_beta    = zeros(N_campioni,1);
T_rec_alpha   = zeros(N_campioni,1);

% Metriche sforzo di controllo
max_u1_cmd = zeros(N_campioni,1);
max_u2_cmd = zeros(N_campioni,1);

max_u1_sat = zeros(N_campioni,1);
max_u2_sat = zeros(N_campioni,1);

%% 3. FIGURE

figure('Name','Monte Carlo Simulink LQG');

subplot(3,1,1)
hold on
grid on
title('Tracking \alpha')
xlabel('Tempo [s]')
ylabel('\alpha [rad]')

subplot(3,1,2)
hold on
grid on
title('Tracking \beta')
xlabel('Tempo [s]')
ylabel('\beta [rad]')

subplot(3,1,3)
hold on
grid on
title('Sforzo di controllo')
xlabel('Tempo [s]')
ylabel('u')

%% ==========================================
% METRICHE
% ==========================================

% Istante del disturbo
t_dist = 10;

% Soglie di recupero
soglia_beta  = 0.07;   % rad
soglia_alpha = 0.01;   % rad

finestra_rec = 1.0;   % deve rimanere nella fascia per 1 s

%% ==========================================
% 4. MONTE CARLO
% ==========================================

for i = 1:N_campioni


    %% 4.1 Estrazione campione incerto

    Pk = usample(P_full_unc);

    [Ak,Bk,Ck,Dk] = ssdata(Pk);


    %% 4.2 Passaggio del campione a Simulink

    assignin('base','A_nom',Ak);
    assignin('base','B_nom',Bk);
    assignin('base','Cy',Ck);
    assignin('base','Dy',Dk);


    %% 4.3 Simulazione Simulink

    simOut = sim('sim_elicottero',...
        'SrcWorkspace','base',...
        'ReturnWorkspaceOutputs','on');


    %% 4.4 Estrazione risultati

    alpha_resp = simOut.alpha_LQG_int.Data;
    beta_resp  = simOut.beta_LQG_int.Data;

    t = simOut.alpha_LQG_int.Time;

    u_cmd = simOut.u_cmd.Data;
    u_sat = simOut.u_sat.Data;

    if size(u_cmd,2) ~= 2
        error('u_cmd deve avere 2 colonne: [u1 u2].');
    end
    
    if size(u_sat,2) ~= 2
        error('u_sat deve avere 2 colonne: [u1_sat u2_sat].');
    end


    %% 4.5 Plot

    subplot(3,1,1)
    plot(t,alpha_resp)
    
    subplot(3,1,2)
    plot(t,beta_resp)
    
    subplot(3,1,3)
    plot(t,u_cmd(:,1))
    plot(t,u_cmd(:,2))



    % ==========================================
    % 1. Picco assoluto
    % ==========================================
    
    max_alpha(i) = max(abs(alpha_resp));
    max_beta(i)  = max(abs(beta_resp));
    
    
    % ==========================================
    % 2. RMS rispetto al riferimento
    % ==========================================
    
    rms_alpha(i) = rms(alpha_resp - 0.1);
    rms_beta(i)  = rms(beta_resp);
    
    
    % ==========================================
    % 3. Errore finale
    % Media degli ultimi 2 s per attenuare il rumore
    % ==========================================
    
    idx_final = t >= (t(end)-2);
    
    final_alpha(i) = abs(mean(alpha_resp(idx_final)) - 0.1);
    
    beta_mean_ss = mean(beta_resp(idx_final));
    
    final_beta(i) = abs(beta_mean_ss);
    
    beta_dev_ss(i) = max(abs(beta_resp(idx_final) - beta_mean_ss));
        
    
    % ==========================================
    % 4. Tempo recupero beta
    % |beta| <= 0.07 rad per almeno 1 s
    % ==========================================
    
    idx_post = find(t >= t_dist);
    
    beta_post = abs(beta_resp(idx_post));
    
    N_finestra = max(1,round(finestra_rec/mean(diff(t))));
    
    T_rec_beta(i) = NaN;
    
    for k = 1:(length(beta_post)-N_finestra+1)
    
        finestra = beta_post(k:k+N_finestra-1);
    
        if all(finestra <= soglia_beta)
    
            T_rec_beta(i) = t(idx_post(k)) - t_dist;
            break
    
        end
    end
    
    if isnan(T_rec_beta(i))
        T_rec_beta(i) = t(end)-t_dist;
    end
    
    
    % ==========================================
    % 5. Tempo recupero alpha
    % |alpha-0.1| <= 0.01 rad per almeno 1 s
    % ==========================================
    
    errore_alpha_post = abs(alpha_resp(idx_post) - 0.1);
    
    T_rec_alpha(i) = NaN;
    
    for k = 1:(length(errore_alpha_post)-N_finestra+1)
    
        finestra = errore_alpha_post(k:k+N_finestra-1);
    
        if all(finestra <= soglia_alpha)
    
            T_rec_alpha(i) = t(idx_post(k)) - t_dist;
            break
    
        end
    end
    
    if isnan(T_rec_alpha(i))
        T_rec_alpha(i) = t(end)-t_dist;
    end

    %% ==========================================
    % 6. SFORZO DI CONTROLLO
    % ==========================================
    
    max_u1_cmd(i) = max(abs(u_cmd(:,1)));
    max_u2_cmd(i) = max(abs(u_cmd(:,2)));
    
    max_u1_sat(i) = max(abs(u_sat(:,1)));
    max_u2_sat(i) = max(abs(u_sat(:,2)));

    % ==========================================
    % Stampa campione
    % ==========================================
    
    fprintf('Campione %d/%d\n',i,N_campioni);
    fprintf('   max |alpha| = %.6f rad\n',max_alpha(i));
    fprintf('   max |beta|  = %.6f rad\n',max_beta(i));
    fprintf('   err finale alpha = %.6f rad\n',final_alpha(i));
    fprintf('   err finale beta  = %.6f rad\n',final_beta(i));
    fprintf('   T recupero beta  = %.3f s\n',T_rec_beta(i));
    fprintf('   T recupero alpha = %.3f s\n',T_rec_alpha(i));
    fprintf('   max |u1_cmd| = %.6f\n',max_u1_cmd(i));
    fprintf('   max |u2_cmd| = %.6f\n',max_u2_cmd(i));
end


%% ==========================================
% 5. RIPRISTINO MODELLO NOMINALE
% ==========================================

assignin('base','A_nom',A_nom_saved);
assignin('base','B_nom',B_nom_saved);
assignin('base','Cy',Cy_saved);
assignin('base','Dy',Dy_saved);


%% ==========================================
% 6. SIMULAZIONE NOMINALE
% ==========================================

disp(' ')
disp('--- Simulazione nominale ---')

simOut_nom = sim('sim_elicottero',...
    'SrcWorkspace','base',...
    'ReturnWorkspaceOutputs','on');


t_nom     = simOut_nom.alpha_LQG_int.Time;
alpha_nom = simOut_nom.alpha_LQG_int.Data;
beta_nom  = simOut_nom.beta_LQG_int.Data;


u_cmd_nom = simOut_nom.u_cmd.Data;
u_sat_nom = simOut_nom.u_sat.Data;

if size(u_cmd_nom,2) ~= 2
    error('u_cmd_nom deve avere 2 colonne: [u1 u2].');
end

if size(u_sat_nom,2) ~= 2
    error('u_sat_nom deve avere 2 colonne: [u1_sat u2_sat].');
end




%% 7. PLOT NOMINALE

subplot(3,1,1)
plot(t_nom,alpha_nom,'k','LineWidth',2)
legend('Campioni incerti','Nominale','Location','best')

subplot(3,1,2)
plot(t_nom,beta_nom,'k','LineWidth',2)
legend('Campioni incerti','Nominale','Location','best')

subplot(3,1,3)
plot(t_nom,u_cmd_nom(:,1),'k','LineWidth',2)
plot(t_nom,u_cmd_nom(:,2),'--k','LineWidth',2)

legend('u_1 campioni','u_2 campioni',...
       'u_1 nominale','u_2 nominale',...
       'Location','best')

fprintf('\n============================================\n')
fprintf('       MONTE CARLO SIMULINK - %d CAMPIONI\n',N_campioni)
fprintf('============================================\n')

fprintf('\n--- ALPHA ---\n')

fprintf('Picco massimo |alpha|       = %.6f rad (campione %d)\n',...
    max(max_alpha),find(max_alpha == max(max_alpha),1));

fprintf('Errore finale massimo alpha = %.6f rad (campione %d)\n',...
    max(final_alpha),find(final_alpha == max(final_alpha),1));

fprintf('RMS massimo errore alpha    = %.6f rad (campione %d)\n',...
    max(rms_alpha),find(rms_alpha == max(rms_alpha),1));

fprintf('Tempo recupero alpha worst  = %.3f s (campione %d)\n',...
    max(T_rec_alpha),find(T_rec_alpha == max(T_rec_alpha),1));


fprintf('\n--- BETA ---\n')

fprintf('Picco massimo |beta|        = %.6f rad (campione %d)\n',...
    max(max_beta),find(max_beta == max(max_beta),1));

fprintf('Errore finale massimo beta  = %.6f rad (campione %d)\n',...
    max(final_beta),find(final_beta == max(final_beta),1));

fprintf('Oscillazione residua beta   = %.6f rad (campione %d)\n',...
    max(beta_dev_ss),find(beta_dev_ss == max(beta_dev_ss),1));

fprintf('RMS massimo beta            = %.6f rad (campione %d)\n',...
    max(rms_beta),find(rms_beta == max(rms_beta),1));

fprintf('Tempo recupero beta worst   = %.3f s (campione %d)\n',...
    max(T_rec_beta),find(T_rec_beta == max(T_rec_beta),1));


fprintf('\n--- SFORZO DI CONTROLLO ---\n')

fprintf('Worst max |u1_cmd| = %.6f (campione %d)\n',...
    max(max_u1_cmd),find(max_u1_cmd == max(max_u1_cmd),1));

fprintf('Worst max |u2_cmd| = %.6f (campione %d)\n',...
    max(max_u2_cmd),find(max_u2_cmd == max(max_u2_cmd),1));

fprintf('Worst max |u1_sat| = %.6f (campione %d)\n',...
    max(max_u1_sat),find(max_u1_sat == max(max_u1_sat),1));

fprintf('Worst max |u2_sat| = %.6f (campione %d)\n',...
    max(max_u2_sat),find(max_u2_sat == max(max_u2_sat),1));

fprintf('\nUtilizzo massimo attuatore:\n')
fprintf('u1: %.2f %% del limite\n',...
    100*max(max_u1_cmd)/umax);

fprintf('u2: %.2f %% del limite\n',...
    100*max(max_u2_cmd)/umax);

n_sat_u1 = sum(max_u1_cmd > umax);
n_sat_u2 = sum(max_u2_cmd > umax);

fprintf('\nControllo effettiva saturazione:\n');

fprintf('max |u1_cmd| - max |u1_sat| = %.6f\n',...
    max(max_u1_cmd) - max(max_u1_sat));

fprintf('max |u2_cmd| - max |u2_sat| = %.6f\n',...
    max(max_u2_cmd) - max(max_u2_sat));

fprintf('\nCampioni con saturazione u1: %d/%d\n',...
    n_sat_u1,N_campioni);

fprintf('Campioni con saturazione u2: %d/%d\n',...
    n_sat_u2,N_campioni);

%% ==========================================
% DISTRIBUZIONE DELLE METRICHE MONTE CARLO
% ==========================================

figure('Name','Distribuzione metriche Monte Carlo');

subplot(2,2,1)
histogram(max_beta,10)
grid on
xlabel('max |\beta| [rad]')
ylabel('Numero campioni')
title('Distribuzione picco |\beta|')

subplot(2,2,2)
histogram(T_rec_beta,10)
grid on
xlabel('Tempo recupero [s]')
ylabel('Numero campioni')
title('Distribuzione T_{rec} \beta')

subplot(2,2,3)
histogram(final_beta,10)
grid on
xlabel('Errore medio finale |\beta| [rad]')
ylabel('Numero campioni')
title('Distribuzione errore finale \beta')

subplot(2,2,4)
histogram(beta_dev_ss,10)
grid on
xlabel('Oscillazione residua [rad]')
ylabel('Numero campioni')
title('Distribuzione oscillazione residua \beta')

figure('Name','Distribuzione metriche Alpha');

subplot(1,3,1)
histogram(max_alpha,10)
grid on
xlabel('max |\alpha| [rad]')
ylabel('Numero campioni')
title('Picco |\alpha|')

subplot(1,3,2)
histogram(T_rec_alpha,10)
grid on
xlabel('Tempo recupero [s]')
ylabel('Numero campioni')
title('T_{rec} \alpha')

subplot(1,3,3)
histogram(final_alpha,10)
grid on
xlabel('Errore finale |\alpha| [rad]')
ylabel('Numero campioni')
title('Errore finale \alpha')

figure('Name','Distribuzione sforzo di controllo');

subplot(2,2,1)
histogram(max_u1_cmd,10)
grid on
xlabel('max |u_1|')
ylabel('Numero campioni')
title('Picco comando u_1')

subplot(2,2,2)
histogram(max_u2_cmd,10)
grid on
xlabel('max |u_2|')
ylabel('Numero campioni')
title('Picco comando u_2')

subplot(2,2,3)
histogram(max_u1_sat,10)
grid on
xlabel('max |u_{1,sat}|')
ylabel('Numero campioni')
title('Picco u_1 dopo saturazione')

subplot(2,2,4)
histogram(max_u2_sat,10)
grid on
xlabel('max |u_{2,sat}|')
ylabel('Numero campioni')
title('Picco u_2 dopo saturazione')