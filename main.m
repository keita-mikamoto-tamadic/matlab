clear;
close;

% ==============================
% 基本設定
% ==============================
iflinear = true;
control_sw = true;


% syms t theta(t) phi(t) tau(t) g
% % 振り子
% syms Mp Jp l u_1
% % タイヤ
% syms Mw Jw a u_2

% パラメータ
g = 9.81;

Mp = 0.263;
l = 0.426;
Jp = 7.669e-4;
u_1 = 0;
u_3= 0;


Mw = 0.028;
a = 0.058;
Jw = 1.133e-5;
u_2 = 0;
u_4 = 0;


% 運動方程式を立てる
run('createEOM.m');

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
x0 = [0.3, 0, 0.03, 0];


% データ長
n_data = length(t);

% 状態
X = nan(4, n_data);
X(:, 1) = x0;

% 入力
U = nan(2, n_data);
F = nan(2, n_data);

% 制御用パラメータ
x_ref = 0;
xki = nan(1, n_data+1);



kp = 15;
ki = 0.3;

if control_sw == false
    kp = 0;
    ki = 0;
end


% シミュレーション回す
for i_data = 1 : n_data

    % 外力
    if i_data < 1000
        d = 0;
    elseif i_data < 1200
        d = -1;
    else
        d = 0;
    end

    U(:, i_data) = [0, 0];
    

    % X(1) = theta, X(2) = dtheta, X(3) = phi, X(4) = dphi
    % X(:, i_data+1) = f_rk4(@state_eq, X(:,i_data), -U(:,i_data), F(:, i_data), delT);
    X(:, i_data+1) = f_rk4(@state_eq, X(:,i_data), U(:,i_data), delT);

end

% 求めた角度から振子位置を計算
X_dec = nan(5, n_data+1);
for i_data = 1:n_data+1

    % theta
    X_dec(:, i_data) = x2rad(X(1,i_data), X(3,i_data), a, l);
end


% ==============================
% アニメーション
% ==============================
hight =400;
width =400;

sizen = 512;



FileName = 'animation.gif';
interval = 30;
hFig = figure('Position',[100 100 width hight]);
htext = [];

axis tight manual

limit = [-1,1];
xlim([limit(1), limit(2)]);
ylim([limit(1), limit(2)]);

axis square;

hold on;

lines = gobjects(4,1);
lines(1) = plot(nan, nan, 'o', 'MarkerSize', 15);
lines(2) = plot(nan, nan, 'o', 'MarkerSize', 1);
lines(3) = plot(nan, nan, '-', 'LineWidth', 1);
lines(4) = plot(nan, nan, '-', 'LineWidth', 1);
 
plot([-5, 5], [-a, -a], '-', 'Color', 'black');

for i_data = 1:interval:n_data+1

    delete(htext);

    x_p = X_dec(2, i_data);
    y_p = X_dec(3, i_data);
    x_w = X_dec(1, i_data);
    x_phi = X_dec(4, i_data);
    y_phi = X_dec(5, i_data);

    set(lines(1), 'XData', x_w, 'YData', 0);
    set(lines(2),'XData', x_p,'YData', y_p);
    set(lines(3),'XData', [x_p, x_w],'YData', [y_p, 0]);
    set(lines(4),'XData', [x_phi, x_w],'YData', [y_phi, 0]);
    htext = text(-0.9, 0.9, ['Time: ', num2str(i_data*delT, 2),'[s]'], 'FontSize', 12);
    
    drawnow;

    % GIFファイルにフレームを追加
    frameData = getframe(hFig);
    im = frame2im(frameData);
    [imind,cm] = rgb2ind(im, sizen);

    % GIFファイルにフレームを追加
    if i_data == 1
         imwrite(imind,cm,FileName,'gif', 'Loopcount',inf,'DelayTime', delT*interval);
    else
        imwrite(imind,cm,FileName,'gif','WriteMode','append','DelayTime',delT*interval);
    end
    pause(delT*interval*0.01);
end
