%% LQG Servo Controller: 1-DOF, 2-DOF, and h2syn equivalents
%
% Plant: G(s) = 3*(-2s+1) / ((5s+1)(10s+1))
%
% Three equivalent synthesis approaches compared:
%   1) Manual LQG: LQR on augmented plant + Kalman filter on original plant
%   2) lqg() command: MATLAB built-in with '1dof'/'2dof' option
%   3) h2syn(): H2 synthesis via generalized plant Pss
%
% The three approaches should give identical closed-loop responses.
%
% Controller structures:
%   1-DOF: K receives e = r-y (tracking error) + xI = integral(e)
%   2-DOF: K receives [e; y; xI] separately (feedforward + feedback)
%
% Key design choice: weight only on xI (integral state), not on x.
% This guarantees zero steady-state error for constant references.

close all; 
clear all; 
clc;

%% =========================================================
%  PLANT MODEL AND DIMENSIONS
%% =========================================================
G = tf(3*[-2 1], conv([5 1],[10 1]));   % G(s) = 3(-2s+1)/((5s+1)(10s+1))
[a,b,c,d] = ssdata(G);
Plant = ss(a,b,c,d);

p      = size(c,1);      % number of outputs
[n, m] = size(b);        % number of states and inputs

% Convenience zero matrices
Znm = zeros(n,m);  Zmm = zeros(m,m);
Znn = zeros(n,n);  Zmn = zeros(m,n);

fprintf('Plant order: %d, inputs: %d, outputs: %d\n', n, m, p);

%% =========================================================
%  STEP 1: LQR STATE FEEDBACK ON AUGMENTED PLANT
%
%  Augmented state: xa = [x(n); xI(m)]
%    xI_dot = r - y = -c*x   (tracking error integrator, r=0 in design)
%
%  Augmented dynamics:
%    xa_dot = Aa*xa + Ba*u
%    Aa = [a,  0 ]    Ba = [b ]
%         [-c, 0 ]         [-d]
%
%  LQR cost: penalize xI only (not x) to ensure integral action
%    J = integral(xa'*Q*xa + u'*R*u) dt
%    Q = blkdiag(0, I)  ->  only xI penalized
%
%  Result: Kra = [Krp | Kri]
%    u* = -Krp*x - Kri*xI   (state feedback + integral action)
%% =========================================================
Aa = [a,  Znm; -c, Zmm];       % augmented plant dynamics
Ba = [b;  -d];                 % augmented plant input matrix

Q  = blkdiag(Znn, eye(m));     % cost: weight on xI only
R  = eye(m);                   % input penalty

Kra = lqr(Aa, Ba, Q, R);
Krp = Kra(1:m, 1:n);           % state feedback gain   (m x n)
Kri = Kra(1:m, n+1:n+m);       % integral feedback gain (m x m)

fprintf('State feedback gain Krp =\n'); disp(Krp);
fprintf('Integral gain      Kri =\n'); disp(Kri);

%% =========================================================
%  STEP 2: KALMAN FILTER ON ORIGINAL PLANT
%
%  The Kalman filter estimates x only (NOT xI).
%  Reason: xI = integral(r-y) is computed directly from measurements
%  and is not a state of the physical plant. Estimating it via Kalman
%  would require (Aa, c_aug) to be detectable, but the s=0 mode of
%  xI is unobservable from y = c*x (the c_aug row for xI is zero).
%
%  Noise model options (switch comment to select):
%    a) Process noise enters all states: Bnoise = eye(n)
%    b) Actuation noise enters through b: Bnoise = b
%% =========================================================

% --- Process noise model (default) ---
Bnoise = eye(n);               % noise enters all states equally
W      = eye(size(Bnoise,2));  % process noise covariance
V      = eye(p);               % measurement noise covariance

% % --- Actuation noise model (uncomment to use) ---
% Bnoise = b;
% W      = eye(m);
% V      = eye(p);

% Estimation model: inputs [u; w_proc], output y
Estss = ss(a, [b, Bnoise], c, zeros(p, m+size(Bnoise,2)));
[~, Ke, ~] = kalman(Estss, W, V);   % Ke: optimal Kalman gain (n x p)

fprintf('Kalman gain Ke =\n'); disp(Ke);

%% =========================================================
%  STEP 3: CONTROLLER REALIZATION
%
%  Controller states: xc = [xI(m); xhat(n)]
%  (integrator first, then observer states)
%
%  Dynamics:
%    xI_dot   = r - y
%    xhat_dot = (a - b*Krp - Ke*c)*xhat - b*Kri*xI + Ke*y
%
%  Control law:
%    u = -Krp*xhat - Kri*xI
%
%  State-space matrices:
%    Ac  = [    0,          0    ]   (m+n) x (m+n)
%          [-b*Kri,  a-b*Krp-Ke*c]
%
%    Bcr = [I ]  input from r   (m+n) x m
%          [0 ]
%
%    Bcy = [-I]  input from y   (m+n) x m
%          [Ke]
%
%    Cc  = [-Kri, -Krp]          m x (m+n)
%    Dc  = 0
%% =========================================================
Ac  = [Zmm,    Zmn;
      -b*Kri,  a - b*Krp - Ke*c];

Bcr = [eye(m); Znm];    % r -> xI_dot = +r,  xhat_dot = 0
Bcy = [-eye(m); Ke ];   % y -> xI_dot = -y,  xhat_dot = Ke*y

Cc  = [-Kri, -Krp];     % u = -Kri*xI - Krp*xhat
Dc  = zeros(m, m);

% --- 1-DOF: K receives -y, r enters through feedback loop ---
% Sign convention: input is -y so that u = K*(-y) = -K*y.
% When loop is closed (y = G*u + r), the integrator sees r-y implicitly.
Klqg = ss(Ac, -Bcy, Cc, -Dc, ...
          'InputName', {'yneg'}, 'OutputName', {'u'});

% --- 2-DOF: K receives [r; y] separately ---
% r feeds the integrator directly (feedforward path).
% y feeds the Kalman observer and the integrator (feedback path).
Klqg2 = ss(Ac, [Bcr, Bcy], Cc, [Dc, Dc], ...
           'InputName', {'r','y'}, 'OutputName', {'u'});

fprintf('Manual controller order: %d\n', size(Ac,1));

%% =========================================================
%  CLOSED-LOOP SIMULATION — MANUAL DESIGN
%% =========================================================

% 1-DOF: standard negative feedback  r -> (+) -> [G*K] -> y -> (-)
CL1  = feedback(G * Klqg, 1);         % r -> y
CL1u = feedback(Klqg, G);             % r -> u

% 2-DOF: sign of feedback already in Bcy, use +1 in feedback()
% feedback(sys, K, feedin, feedout, sign):
%   feedin=2 (y input of K), feedout=1 (y output of G), sign=+1
CL2   = feedback(G * Klqg2, 1, 2, 1, +1);
sys2  = CL2 * [1; 0];    % r -> y  (select r channel)
CL2u  = feedback(Klqg2, G, 2, 1, +1);
sys2u = CL2u * [1; 0];   % r -> u

figure('Name','LQG Servo: Manual Design','Position',[50 50 900 600]);
subplot(2,1,1);
step(CL1); hold on; step(sys2);
yline(1,'k--','LineWidth',1);
legend('1-DOF','2-DOF','r=1');
ylabel('y(t)'); title('Output y'); grid on;
subplot(2,1,2);
step(CL1u); hold on; step(sys2u);
yline(0,'k--','LineWidth',1);
legend('1-DOF','2-DOF');
ylabel('u(t)'); title('Control effort u'); grid on;
sgtitle('LQG Servo: Manual Design','FontSize',13);

%% =========================================================
%  ALTERNATIVE: lqg() COMMAND
%
%  lqg() builds the LQG servo controller internally.
%  QXU = blkdiag(Qx, R): weights on [x; u] for the original plant
%  QWV = blkdiag(W, V):  noise covariances
%  QI:                   weight on integral states (= eye(m) here)
%% =========================================================
QXU = blkdiag(Q(1:n,1:n), R);    % [x; u] weights (Qx=0, R=I)
QWV = blkdiag(W, V);              % noise covariances
QI  = Q(n+1:end, n+1:end);       % integral weight = eye(m)

[reg2, ~] = lqg(Plant, QXU, QWV, QI, '2dof');
[reg1, ~] = lqg(Plant, QXU, QWV, QI, '1dof');

fprintf('\nlqg() 1-DOF order: %d,  2-DOF order: %d\n', order(reg1), order(reg2));

CL1_lqg  = feedback(G * reg1, 1);
CL2_lqg  = feedback(G * reg2, 1, 2, 1, +1);
sys2_lqg = CL2_lqg * [1; 0];

figure('Name','Manual vs lqg()','Position',[50 50 900 500]);
step(CL1); hold on; step(CL1_lqg,'--');
step(sys2);         step(sys2_lqg,'--');
legend('1-DOF (manual)','1-DOF (lqg)','2-DOF (manual)','2-DOF (lqg)');
title('Manual LQG vs lqg() — should overlap'); grid on;

%% =========================================================
%  ALTERNATIVE: H2SYN, i.e., LQG VIA GENERALIZED PLANT
%
%  h2syn minimizes ||lft(Pss,K)||_2, equivalent to LQG when
%  Pss correctly encodes Q, R, W, V.
%
%  Correspondence LQG <-> H2:
%    Regulated outputs: z = [Q^{1/2}*xa; R^{1/2}*u]
%    Exogenous inputs:  w = [W^{1/2}*w_proc; V^{1/2}*w_meas]
%
%  KEY CHALLENGE: the integrator state xI has a pole at s=0.
%  h2syn requires (Aa, C2) to have no unobservable modes on jw-axis.
%  Since C2 = [c, 0] (xI column is zero), xI is NOT observable from y.
%
%  SOLUTION: add xI as an explicit measured output.
%  xI is the internal integral state — it is available to the controller.
%  With xI measured directly, the s=0 mode becomes observable and
%  h2syn can find the optimal controller with integral action.
%
%  INPUT ORDERING (required by h2syn): [...exogenous w..., u]  (u LAST)
%  OUTPUT ORDERING (required by h2syn): [z..., v...]  (z first, v last)
%% =========================================================
nw = size(Bnoise, 2);   % number of process noise channels

% Augmented plant matrices (pure integrator, no pole shift needed)
Aa_h2 = [a,  Znm; -c, Zmm];         % xI_dot = -c*x = -y  (r=0 for 1DOF)
Ba_u  = [b;  -d  ];                  % control input channel
Ba_r  = [zeros(n,m); eye(m)];        % reference input channel (xI_dot += r)
Ba_w  = [Bnoise; zeros(m,nw)];       % process noise channel (on x only)

% Output weighting matrices
Qsqrt = blkdiag(zeros(n), eye(m));  % Q^{1/2}: penalize xI only
Rsqrt = eye(m);                      % R^{1/2}: penalize u

% Plant output and integral state observation matrices
Ca_y  = [c, zeros(p,m)];            % y  = c*x        (no xI term)
Ca_xI = [zeros(m,n), eye(m)];       % xI = last m states of xa

%% --- 2-DOF Pss ---
%
%  K receives [e; y; xI] where:
%    e  = r - y + w_meas  (tracking error, measured with noise)
%    y  = c*x  + w_meas   (plant output, measured with noise)
%    xI = integral(r-y)   (integral state, directly observable)
%
%  All three measured outputs ensure (Aa_h2, C2) has no jw-axis
%  unobservable modes, satisfying h2syn conditions.
%
%  Inputs:  [r(m), w_proc(nw), w_meas(p), u(m)]
%  Outputs: [y_out(p), z_xa(n+m), z_u(m), e(p), y(p), xI(m)]
%    y_out: added as first regulated output to survive lft() and
%           allow extraction of r->y transfer function.

ni2 = m + nw + p + m;

B2 = [Ba_r, Ba_w, zeros(n+m,p), Ba_u];

C2 = [Ca_y;               % y_out  (regulated, preserved after lft)
      Qsqrt;               % z_xa = Q^{1/2}*xa
      zeros(m,n+m);        % z_u  (feedthrough via D only)
      -Ca_y;               % e    = -c*x + r + w_meas  (= r-y+noise)
      Ca_y;                % y    =  c*x     + w_meas
      Ca_xI];              % xI   (no noise, direct state readout)

%     r(m)          w_proc(nw)        w_meas(p)      u(m)
D2 = [zeros(p,m),   zeros(p,nw),     zeros(p,p),    zeros(p,m);    % y_out
      zeros(n+m,m), zeros(n+m,nw),   zeros(n+m,p),  zeros(n+m,m); % z_xa
      zeros(m,m),   zeros(m,nw),     zeros(m,p),    Rsqrt;          % z_u
      eye(m),       zeros(p,nw),     eye(p),         zeros(p,m);    % e = r-y+noise
      zeros(p,m),   zeros(p,nw),     eye(p),         zeros(p,m);    % y = noise
      zeros(m,m),   zeros(m,nw),     zeros(m,p),    zeros(m,m)];   % xI

Pss2_ext = ss(Aa_h2, B2, C2, D2);
nmeas2   = p + p + m;    % K measures [e; y; xI]

fprintf('\n--- h2syn 2-DOF ---\n');
[Kh2_2, ~, gam_h2_2] = h2syn(Pss2_ext, nmeas2, m);
fprintf('H2 norm: %.4f,  controller order: %d\n', gam_h2_2, order(Kh2_2));

% lft(Pss2_ext, Kh2_2) closes the [e;y;xI] <-> u loop.
% Remaining inputs:  [r, w_proc, w_meas]
% Remaining outputs: [y_out, z_xa, z_u]
CL_h2_2 = lft(Pss2_ext, Kh2_2);
sys_h2_2 = CL_h2_2(1, 1);    % channel r -> y_out
fprintf('2-DOF h2syn DC gain = %.4f  (should be 1)\n', dcgain(sys_h2_2));

%% --- 1-DOF Pss ---
%
%  K receives [e; xI] where:
%    e  = r - y + w_meas  (tracking error)
%    xI = integral(r-y)   (integral state, directly observable)
%
%  Note: compared to 2-DOF, K does NOT receive y separately.
%  This is the 1-DOF structure: only the error (and its integral)
%  are available to the controller, not r and y individually.
%
%  Inputs:  [r(m), w_proc(nw), w_meas(p), u(m)]
%  Outputs: [y_out(p), z_xa(n+m), z_u(m), e(p), xI(m)]

B1e = [Ba_r, Ba_w, zeros(n+m,p), Ba_u];

C1e = [Ca_y;               % y_out  (regulated, preserved after lft)
       Qsqrt;               % z_xa
       zeros(m,n+m);        % z_u (via D)
       -Ca_y;               % e  = -c*x + r + w_meas
       Ca_xI];              % xI (direct state readout)

%      r(m)          w_proc(nw)        w_meas(p)      u(m)
D1e = [zeros(p,m),   zeros(p,nw),     zeros(p,p),    zeros(p,m);    % y_out
       zeros(n+m,m), zeros(n+m,nw),   zeros(n+m,p),  zeros(n+m,m); % z_xa
       zeros(m,m),   zeros(m,nw),     zeros(m,p),    Rsqrt;          % z_u
       eye(p),       zeros(p,nw),     eye(p),         zeros(p,m);    % e = r-y+noise
       zeros(m,m),   zeros(m,nw),     zeros(m,p),    zeros(m,m)];   % xI

Pss1e   = ss(Aa_h2, B1e, C1e, D1e);
nmeas1e = p + m;    % K measures [e; xI]

fprintf('\n--- h2syn 1-DOF ---\n');
[Kh2_1, ~, gam_h2_1] = h2syn(Pss1e, nmeas1e, m);
fprintf('H2 norm: %.4f,  controller order: %d\n', gam_h2_1, order(Kh2_1));

% lft(Pss1e, Kh2_1) closes [e;xI] <-> u loop.
% Remaining inputs:  [r, w_proc, w_meas]
% Remaining outputs: [y_out, z_xa, z_u]
CL_h2_1 = lft(Pss1e, Kh2_1);
sys_h2_1 = CL_h2_1(1, 1);    % channel r -> y_out
fprintf('1-DOF h2syn DC gain = %.4f  (should be 1)\n', dcgain(sys_h2_1));

%% =========================================================
%  COMPARISON PLOT: Manual LQG vs lqg() vs h2syn
%
%  All three designs use identical Q, R, W, V weights.
%  Curves should overlap if the three formulations are equivalent.
%% =========================================================
figure('Name','Comparison: Manual vs lqg() vs h2syn','Position',[50 50 950 550]);
t_cmp = 0:0.01:60;

plot(t_cmp, step(sys2,     t_cmp), 'b-',  'LineWidth',1.5); hold on;
plot(t_cmp, step(sys2_lqg, t_cmp), 'r--', 'LineWidth',1.5);
plot(t_cmp, step(CL1,      t_cmp), 'c-.',  'LineWidth',1.5);
plot(t_cmp, step(CL1_lqg,  t_cmp), 'k--', 'LineWidth',1.5);
leg = {'2-DOF manual','2-DOF lqg()','1-DOF manual','1-DOF lqg()'};

if exist('sys_h2_2','var')
    plot(t_cmp, step(sys_h2_2, t_cmp), 'g:', 'LineWidth',2);
    leg{end+1} = '2-DOF h2syn';
end
if exist('sys_h2_1','var')
    plot(t_cmp, step(sys_h2_1, t_cmp), 'm:', 'LineWidth',2);
    leg{end+1} = '1-DOF h2syn';
end

yline(1,'k--','LineWidth',1,'HandleVisibility','off');
legend(leg,'Location','best');
title('Manual LQG vs lqg() vs h2syn (should overlap)');
xlabel('Time (s)'); ylabel('y(t)'); grid on;