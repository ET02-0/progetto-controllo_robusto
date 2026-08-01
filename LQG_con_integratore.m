% =========================================================================
% Script di Inizializzazione e Linearizzazione: Elicottero 2DoF
% Progetto di Controllo Robusto (Sintesi Robusta su 4 stati)
% =========================================================================

close all
clc

disp('SINTESI LQG 2DOF CON INTEGRATORE')


%% Plant nominale

[A_heli,B_heli,C_ext,D_ext] = ssdata(P_nom);


%% Attuatori

[A_act,B_act,C_act,D_act] = ssdata(G_act_2nd);

Actuator = ss(A_act,B_act,C_act,D_act);

Actuators_MIMO = blkdiag(Actuator,Actuator);


%% Plant esteso

P_ext = P_nom*Actuators_MIMO;

[A_ext,B_ext,C_ext,D_ext]=ssdata(P_ext);

M = [A_ext B_ext;
     -C_ext zeros(2)];

disp(rank(M))
disp(size(M,1))
tzero(P_ext)

n_ext = size(A_ext,1);
m     = size(B_ext,2);
p     = size(C_ext,1);


disp('Dimensioni modello esteso:')
fprintf('Stati   = %d\n',n_ext)
fprintf('Ingressi= %d\n',m)
fprintf('Uscite  = %d\n',p)


if n_ext ~= 8
    error('Errore: il modello esteso non ha 8 stati.')
end



%% -------------------------------------------------
% 4) Analisi modello esteso
% --------------------------------------------------

disp('Poli modello esteso:')
disp(eig(A_ext))


figure
pzmap(P_ext)
grid on
title('Poli modello esteso')

%% Integratore solo sul pitch

C_int = [1 0 0 0 0 0 0 0];

ni = size(C_int,1);


Aa = [
A_ext zeros(n_ext,ni);
-C_int zeros(ni,ni)
];


Ba = [
B_ext;
zeros(ni,m)
];

rank_ctrb = rank(ctrb(Aa,Ba));

if rank_ctrb < size(Aa,1)

    fprintf("Modo non controllabile:\n")

    [V,D]=eig(Aa);

    for i=1:length(diag(D))
        if abs(real(D(i)))<1e-8
            disp(V(:,i))
        end
    end

end


%% -------------------------------------------------
% 6) Verifica controllabilità
% --------------------------------------------------

rank_ctrb = rank(ctrb(Aa,Ba));

fprintf('\nControllabilita: %d/%d\n',...
    rank_ctrb,size(Aa,1));


if rank_ctrb ~= size(Aa,1)
    warning('Sistema aumentato non completamente controllabile')
end
eig(Aa)
Co = ctrb(Aa,Ba);

[U,S,V]=svd(Co);

disp('Singular values controllabilità:')
disp(diag(S))


%% -------------------------------------------------
% 7) Verifica osservabilità
% --------------------------------------------------

rank_obs = rank(obsv(A_ext,C_ext));

fprintf('Osservabilita: %d/%d\n',...
    rank_obs,size(A_ext,1));


if rank_obs ~= size(A_ext,1)
    warning('Sistema non completamente osservabile')
end



%% -------------------------------------------------
% 8) Pesi LQR
%
% Stato:
% x = [elicottero(4); attuatori(4)]
%
% --------------------------------------------------


Q_heli = diag([
100;   % alpha
600;   % alpha_dot
150;   % beta
250    % beta_dot
]);

Q_act = diag([5,0.5,5,0.5]);

R_lqr = 3*eye(2);


Q_int = 200;


Q_lqr = blkdiag(Q_heli,Q_act,Q_int);



%% -------------------------------------------------
% 9) Calcolo LQR
% --------------------------------------------------

K_aug = lqr(Aa,Ba,Q_lqr,R_lqr);
umax = 5; % esempio Nm

disp('Guadagni LQR')
disp(K_aug)

Acl = Aa-Ba*K_aug;

disp('Poli sistema chiuso LQR:')
disp(eig(Acl))

Krp = K_aug(:,1:n_ext);

Kri = K_aug(:,n_ext+1:end);



disp('Guadagno LQR calcolato')



%% -------------------------------------------------
% 10) Filtro di Kalman
% --------------------------------------------------
Gk = B_ext;


Qn = 5e-3*eye(m);  % rumore processo stati
Rn = 5e-2*eye(p);  % rumore misura


Estimator = ss(...
    A_ext,...
    Gk,...
    C_ext,...
    zeros(p,m));


[~,Ke,~] = kalman(Estimator,Qn,Rn);



disp('Filtro Kalman calcolato')



%% -------------------------------------------------
% 11) Costruzione controllore LQG 2DOF con integratore su alfa
% -------------------------------------------------

ni = size(C_int,1);


Ac_ctrl = [
    zeros(ni,ni)              -C_int;
    -B_ext*Kri                A_ext-B_ext*Krp-Ke*C_ext
];


Bcr = [
    eye(ni);
    zeros(n_ext,ni)
];




Bcy = [
    zeros(ni,p);
    Ke
];


Cc_ctrl = [
    -Kri   -Krp
];


Dc_ctrl = zeros(m,ni+p);



K_lqg_int = ss(...
    Ac_ctrl,...
    [Bcr Bcy],...
    Cc_ctrl,...
    Dc_ctrl);

K_lqg_int.InputName = {
    'r_alpha'
    'alpha'
    'beta'
};

K_lqg_int.OutputName = {
    'u1'
    'u2'
};

disp('==============================================')
disp(' K_lqg_int creato correttamente')
disp(' Inserire il blocco in Simulink')
size(K_lqg_int)
disp('==============================================')

%% ==========================================
% Connessione LQG 2DOF con riferimento
% ==========================================
P_full = P_nom * blkdiag(G_act_nom,G_act_nom);
P = P_full;
K = K_lqg_int;

% nomi segnali

P.InputName  = {'u1','u2'};
P.OutputName = {'alpha','beta'};


K.InputName = {
    'r_alpha'
    'alpha'
    'beta'
};

K.OutputName = {
    'u1'
    'u2'
};

%% Connessione nominale

Pnom = P_ext;

Pnom.InputName = {'u1','u2'};
Pnom.OutputName = {'alpha','beta'};


K.InputName = {
    'r_alpha'
    'alpha'
    'beta'
};

K.OutputName = {
    'u1'
    'u2'
};


CL_nom = connect(Pnom,K,...
    {'r_alpha'},...
    {'alpha','beta'});
dcgain(CL_nom)
figure
step(CL_nom)

grid on

% collegamento feedback
CL = connect(P,K,...
    {'r_alpha'},...
    {'alpha','beta'});

%% ==========================================
% Analisi stabilità
% ==========================================

disp('==============================================')
disp(' ANALISI STABILITA')
disp('==============================================')

% Poli sistema chiuso nominale
poli_CL = pole(CL);

disp('Poli closed-loop:')
disp(poli_CL)


if all(real(poli_CL)<0)
    disp('Sistema nominale stabile')
else
    disp('ATTENZIONE: sistema nominale INSTABILE')
end


% Margine di stabilità classico
figure
pzmap(CL)
grid on
title('Poli Closed Loop LQG')


% Risposta al gradino
figure
step(CL)
grid on
title('Tracking riferimento alpha')


P_unc = P_full_unc;

P_unc.InputName = {
    'u1'
    'u2'
};


P_unc.OutputName = {
    'alpha'
    'beta'
};


CL_unc = connect(P_unc,K,...
    {'r_alpha'},...
    {'alpha','beta'});

%% ==========================================
% Analisi robustezza
% ==========================================

disp('==============================================')
disp(' ANALISI ROBUSTEZZA')
disp('==============================================')

% Controllo nomi ingressi/uscite
CL_unc.InputName
CL_unc.OutputName


% Robust stability margin
stab_margin = robuststab(CL_unc);

disp('Margine robustezza:')
disp(stab_margin)


% Diagramma valori singolari
omega = logspace(-2,3,500);

figure
sigma(CL_unc,omega)
grid on
title('Singular Values Closed Loop incerto')

%% Robust performance

omega = logspace(-2,3,500);

figure
sigma(CL_unc,omega)
grid on
title('Singular Values Closed Loop')

CL_nom = usubs(CL_unc);

figure
step(CL_nom)
grid on
title('Tracking alpha nominale')
%% Analisi prestazioni

S = stepinfo(CL_nom);

disp('Tempo di assestamento alpha/beta:')
disp(S(1).SettlingTime)
disp(S(2).SettlingTime)

disp('Overshoot alpha/beta:')
disp(S(1).Overshoot)
disp(S(2).Overshoot)
figure
step(CL_nom)

grid on
title('Risposta al gradino nominale')
legend('alpha','beta','Location','best')
wcgain(CL_unc)

size(Aa)
size(Ba)
size(K_aug)
rank(ctrb(Aa,Ba))
eig(Aa-Ba*K_aug)