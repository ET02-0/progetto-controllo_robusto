function m = evalLQG(A,B,C,Bd,Ke,Kx,Ki,Hab,Vpsd,dAmp)

n  = size(A,1);
ni = size(Ki,2);
ny = size(C,1);

Hi = Hab(1:ni,:);

% Metriche inizializzate anche in caso di instabilità
m.stable    = false;
m.peakA     = NaN;
m.peakB     = NaN;
m.rmsB      = NaN;
m.rmsU      = NaN;
m.peakU     = NaN;
m.u1Min     = NaN;
m.u1Max     = NaN;
m.u2Peak    = NaN;
m.rmsA  = NaN;
m.rmsU1 = NaN;
m.rmsU2 = NaN;

% Stati: [x; xhat; xi]
Acl = [ A,        -B*Kx,         -B*Ki;
        Ke*C,      A-B*Kx-Ke*C,  -B*Ki;
       -Hi*C,      zeros(ni,n),  zeros(ni) ];

m.stable = all(real(eig(Acl)) < -1e-8);

if ~m.stable
    return;
end

% Disturbi: Bd deve essere n x 2
Bw_d = [Bd;
        zeros(n+ni,size(Bd,2))];

% Rumore IMU: fattorizzazione della covarianza/PSD
Lv = chol(Vpsd,'lower');

Bw_v = [zeros(n,ny);
        Ke;
       -Hi] * Lv;

Bw = [Bw_d, Bw_v];

% Uscite: alpha, beta, comando 1, comando 2
% Si assume x(1)=alpha e x(3)=beta.
Cout = [1, zeros(1,2*n+ni-1);
        0, 0, 1, zeros(1,2*n+ni-3);
        zeros(2,n), -Kx, -Ki];

G = ss(Acl,Bw,Cout,zeros(4,size(Bw,2)));

t = 0:0.005:15;

% Risposta angolare ai due disturbi simultanei
Y1 = reshape(step(dAmp(1)*G(1:2,1),t),numel(t),2);
Y2 = reshape(step(dAmp(2)*G(1:2,2),t),numel(t),2);
Y  = Y1 + Y2;

m.peakA = max(abs(Y(:,1)));
m.peakB = max(abs(Y(:,2)));

% Rumore IMU -> beta e rumore IMU -> comandi
%% Rumore IMU -> angoli e comandi

idxNoise = size(Bd,2) + (1:ny);

m.rmsA  = norm(G(1,idxNoise),2);
m.rmsB  = norm(G(2,idxNoise),2);

m.rmsU1 = norm(G(3,idxNoise),2);
m.rmsU2 = norm(G(4,idxNoise),2);

m.rmsU = hypot(m.rmsU1,m.rmsU2);

% Risposta dei comandi ai due disturbi simultanei
U1 = reshape(step(dAmp(1)*G(3:4,1),t),numel(t),2);
U2 = reshape(step(dAmp(2)*G(3:4,2),t),numel(t),2);
U  = U1 + U2;

m.u1Min  = min(U(:,1));
m.u1Max  = max(U(:,1));
m.u2Peak = max(abs(U(:,2)));

% Sforzo normalizzato:
% deltaF1 in [-2.5, 2.5], deltaF2 in [-2.5, 2.5]
% peakU = 1 significa raggiungimento di un limite.
m.peakU = max([0, ...
               m.u1Max/2.5, ...
              -m.u1Min/2.5, ...
               m.u2Peak/2.5]);

end