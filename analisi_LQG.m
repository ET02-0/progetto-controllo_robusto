close all
clc

close all
clc

%% ==========================================
% SELEZIONE CONTROLLORE
% 1 = LQG 1-DOF
% 2 = LQG 2-DOF senza integratore
% 3 = LQG 2-DOF con integratore
%% ==========================================

tipo_LQG = 3;

tipo_test = 2;
%tipo_test = 1 --> risposta libera
%tipo_test = 2 --> tracking
%tipo_test = 3 --> reiezione disturbo e rumore
switch tipo_LQG

    case 1
        
        alpha = out.alpha_LQG1;
        beta  = out.beta_LQG1;

        %% Estrazione dati

        t_alpha = alpha.Time;
        alpha_value = alpha.Data;
        
        t_beta = beta.Time;
        beta_value = beta.Data;

        err_alpha = -alpha_value;
        err_beta  = -beta_value;

        
        nome = 'LQG 1-DOF senza integratore';


    case 2
        
        alpha = out.alpha_LQG2;
        beta  = out.beta_LQG2;
        %% Estrazione dati

        t_alpha = alpha.Time;
        alpha_value = alpha.Data;
        
        t_beta = beta.Time;
        beta_value = beta.Data;

        r_alpha = 0.1;
        r_beta  = 0.1;

        err_alpha = r_alpha - alpha_value;
        err_beta  = r_beta - beta_value;
        
        nome = 'LQG 2-DOF senza integratore';

    case 3
        
        alpha = out.alpha_LQG_int;
        beta  = out.beta_LQG_int;
    
        t_alpha = alpha.Time;
        alpha_value = alpha.Data;
    
        t_beta = beta.Time;
        beta_value = beta.Data;
    
        if tipo_test == 1
            
            % risposta libera
            err_alpha = -alpha_value;
            err_beta  = -beta_value;
    
        elseif tipo_test == 2
            
            % tracking
            r_alpha = 0.1;
            r_beta  = 0;
    
            err_alpha = r_alpha-alpha_value;
            err_beta  = r_beta-beta_value;
    
        elseif tipo_test == 3
            
            % disturbo + rumore
            r_alpha = 0;
            r_beta  = 0;
    
            err_alpha = -alpha_value;
            err_beta  = -beta_value;
    
        end
    
        nome = 'LQG 2-DOF con integratore';
end



%% ==========================================
%  ERRORE MASSIMO E RMS
%% ==========================================


max_error_alpha = max(abs(err_alpha));
max_error_beta  = max(abs(err_beta));


RMS_alpha = rms(err_alpha);
RMS_beta  = rms(err_beta);



%% ==========================================
% TEMPO DI ASSESTAMENTO 2%
% tracking
%% ==========================================
if tipo_test == 1

    % risposta libera
    tol_alpha = 0.02*abs(alpha_value(1));
    tol_beta  = 0.02*abs(beta_value(1));

    signal_alpha = alpha_value;
    signal_beta  = beta_value;


elseif tipo_test == 2

    % tracking riferimento
    tol_alpha = 0.02*max(abs(r_alpha),1e-3);
    tol_beta  = 0.02*max(abs(r_beta),1e-3);

    signal_alpha = err_alpha;
    signal_beta  = err_beta;


elseif tipo_test == 3

    % reiezione disturbo + rumore
    tol_alpha = NaN;
    tol_beta  = NaN;

    signal_alpha = err_alpha;
    signal_beta  = err_beta;

end

Ts_alpha = NaN;

for k = 1:length(signal_alpha)

    if all(abs(signal_alpha(k:end)) <= tol_alpha)

        Ts_alpha = t_alpha(k);
        break

    end

end


Ts_beta = NaN;

for k = 1:length(signal_beta)

    if all(abs(signal_beta(k:end)) <= tol_beta)

        Ts_beta = t_beta(k);
        break

    end

end



%% ==========================================
% OVERSHOOT
%% ==========================================

if tipo_test == 1

    alpha_peak = max(abs(alpha_value(2:end)));
    beta_peak  = max(abs(beta_value(2:end)));
    
    overshoot_alpha = max(0,(alpha_peak-abs(alpha_value(1)))/abs(alpha_value(1))*100);
    overshoot_beta  = max(0,(beta_peak-abs(beta_value(1)))/abs(beta_value(1))*100);


elseif tipo_test == 2


    if abs(r_alpha)>1e-6
        overshoot_alpha = (max(alpha_value)-r_alpha)/abs(r_alpha)*100;
    else
        overshoot_alpha = NaN;
    end


    if abs(r_beta)>1e-6
        overshoot_beta = (max(beta_value)-r_beta)/abs(r_beta)*100;
    else
        overshoot_beta = NaN;
    end


elseif tipo_test == 3

    % disturbo + rumore
    overshoot_alpha = NaN;
    overshoot_beta  = NaN;

end


%% ==========================================
%  STEADY STATE CON DISTURBO + RUMORE
%
%  Non usare il 2% perché senza integratore
%  rimane un offset statico
%% ==========================================


N = length(beta_value);


alpha_ss = alpha_value(round(0.8*N):end);
beta_ss  = beta_value(round(0.8*N):end);



mean_alpha_ss = mean(alpha_ss);
mean_beta_ss  = mean(beta_ss);


std_alpha_ss = std(alpha_ss);
std_beta_ss  = std(beta_ss);



max_alpha_ss = max(abs(alpha_ss));
max_beta_ss  = max(abs(beta_ss));



%% ==========================================
%  RISULTATI
%% ==========================================

fprintf('\n=====================================\n')
fprintf('       %s\n',nome)
fprintf('=====================================\n\n')


fprintf('--- Transitorio ---\n')

fprintf('Errore massimo alpha = %.5f rad\n',max_error_alpha);
fprintf('Errore massimo beta  = %.5f rad\n',max_error_beta);

fprintf('RMS alpha = %.5f rad\n',RMS_alpha);
fprintf('RMS beta  = %.5f rad\n',RMS_beta);


fprintf('\nTempo assestamento alpha = %.3f s\n',Ts_alpha);
fprintf('Tempo assestamento beta  = %.3f s\n',Ts_beta);


fprintf('\nOvershoot alpha = %.2f %%\n',overshoot_alpha);
fprintf('Overshoot beta  = %.2f %%\n',overshoot_beta);


if tipo_test == 3
    fprintf('\n--- Reiezione disturbo + rumore ---\n')
    
    
    fprintf('Media alpha steady state = %.5f rad\n',mean_alpha_ss);
    fprintf('Media beta steady state  = %.5f rad\n',mean_beta_ss);
    
    
    fprintf('Deviazione alpha steady state = %.5f rad\n',std_alpha_ss);
    fprintf('Deviazione beta steady state  = %.5f rad\n',std_beta_ss);
    
    
    fprintf('Picco alpha steady state = %.5f rad\n',max_alpha_ss);
    fprintf('Picco beta steady state  = %.5f rad\n',max_beta_ss);

end

%% ==========================================
%  GRAFICI
%% ==========================================

figure

subplot(2,1,1)
plot(t_alpha,alpha_value,'LineWidth',1.5)
grid on
xlabel('Tempo [s]')
ylabel('\alpha [rad]')
title([nome ' - Pitch'])


subplot(2,1,2)
plot(t_beta,beta_value,'LineWidth',1.5)
grid on
xlabel('Tempo [s]')
ylabel('\beta [rad]')
title([nome ' - Yaw'])

x_real = out.xout{5}.Values.Data;
x_hat  = out.xout{6}.Values.Data;
x_hat_plant = x_hat(:,1:4);

t = out.xout{6}.Values.Time;

alpha_hat = x_hat_plant(:,1);
beta_hat  = x_hat_plant(:,3);

alpha_real = x_real(:,1);
beta_real  = x_real(:,3);

figure

subplot(2,1,1)
plot(t,alpha_hat)
grid on
ylabel('$\hat{\alpha}$','Interpreter','latex')
title('Stima Kalman $\alpha$','Interpreter','latex')


subplot(2,1,2)
plot(t,beta_hat)
grid on
ylabel('$\hat{\beta}$','Interpreter','latex')
title('Stima Kalman $\beta$','Interpreter','latex')
xlabel('Tempo [s]')

figure

subplot(2,1,1)

alpha_meas = out.alpha_meas_LQG2;
beta_meas  = out.beta_meas_LQG2;

t_meas_alpha = alpha_meas.Time;
t_meas_beta  = beta_meas.Time;

alpha_meas_value = alpha_meas.Data;
beta_meas_value  = beta_meas.Data;

t = out.xout{6}.Values.Time;


subplot(2,1,1)

plot(t_meas_alpha,alpha_meas_value)
hold on
plot(t,alpha_hat,'LineWidth',1.5)

grid on
legend('Misura rumorosa','Stima Kalman')
title('Filtro Kalman - \alpha')


subplot(2,1,2)

plot(t_meas_beta,beta_meas_value)
hold on
plot(t,beta_hat,'LineWidth',1.5)

grid on
legend('Misura rumorosa','Stima Kalman')
title('Filtro Kalman - \beta')

size(K_lqg_2dof.A)


err_K_alpha = alpha_real - alpha_hat;
err_K_beta  = beta_real - beta_hat;


fprintf('RMS stima alpha = %.5f\n',rms(err_K_alpha))
fprintf('RMS stima beta  = %.5f\n',rms(err_K_beta))


t = out.xout{5}.Values.Time;




figure

subplot(2,1,1)

plot(t,alpha_real,'LineWidth',1.5)
hold on
plot(t,alpha_hat,'LineWidth',1.5)

grid on
legend('Stato reale','Stima Kalman')
title('Errore stima alpha')


subplot(2,1,2)

plot(t,beta_real,'LineWidth',1.5)
hold on
plot(t,beta_hat,'LineWidth',1.5)

grid on
legend('Stato reale','Stima Kalman')
title('Errore stima beta')

figure

subplot(2,1,1)

plot(t,err_K_alpha,'LineWidth',1.5)
grid on
xlabel('Tempo [s]')
ylabel('e_\alpha [rad]')
title('Errore filtro Kalman - alpha')


subplot(2,1,2)

plot(t,err_K_beta,'LineWidth',1.5)
grid on
xlabel('Tempo [s]')
ylabel('e_\beta [rad]')
title('Errore filtro Kalman - beta')


std_K_alpha = std(err_K_alpha);
std_K_beta  = std(err_K_beta);

fprintf('STD errore Kalman alpha = %.5f\n',std_K_alpha)
fprintf('STD errore Kalman beta  = %.5f\n',std_K_beta)

max_K_alpha = max(abs(err_K_alpha));
max_K_beta  = max(abs(err_K_beta));

fprintf('MAX errore Kalman alpha = %.5f\n',max_K_alpha)
fprintf('MAX errore Kalman beta  = %.5f\n',max_K_beta)

[max_K_beta,idx] = max(abs(err_K_beta));

fprintf('Picco beta al tempo %.3f s\n',t(idx))

%% RMS Kalman steady state

N = length(err_K_beta);

err_K_alpha_ss = err_K_alpha(round(0.8*N):end);
err_K_beta_ss  = err_K_beta(round(0.8*N):end);


fprintf('\n--- Kalman steady state ---\n')

fprintf('RMS alpha SS = %.5f rad\n',rms(err_K_alpha_ss))
fprintf('RMS beta SS  = %.5f rad\n',rms(err_K_beta_ss))


fprintf('MAX alpha SS = %.5f rad\n',max(abs(err_K_alpha_ss)))
fprintf('MAX beta SS  = %.5f rad\n',max(abs(err_K_beta_ss)))

plot(t,beta_real)
hold on
plot(t,beta_hat)

plot(t,err_K_beta)