clear;
close;

% ==============================
% 基本設定
% ==============================


% syms t theta(t) phi(t) tau(t) g
% % 振り子
% syms Mp Jp l u_1
% % タイヤ
% syms Mw Jw a u_2

% パラメータ
g = 9.81;

Mp = 0.5;
Jp = 0.1;
l = 0.5;
u_1 = 0.24;

Mw = 0.3;
Jw = 1;
a = 0.3;
u_2 = 0.1;

run('LagrangeEom.mlx');

% ==============================
% シミュレーション
% ==============================

% 各種設定値
% 時間設定
delT = 0.001;
t0 = 0;
tf = 10;
t = t0 : delT : tf;

% 初期状態[theta, dtheta, phi, dphi]
x0 = [0.03, 0, 0.03, 0];


% データ長
n_data = length(t);

X = nan(4, n_data);
X(:, 1) = x0;

U = nan(2, n_data);

% シミュレーション回す
for i_data = 1 : n_data


    U(:, i_data) = zeros(2,1);
    X(:, i_data+1) = f_rk4(@state_eq, X(:,i_data), U(:,i_data), delT);
    
end