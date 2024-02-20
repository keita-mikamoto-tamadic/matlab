function x_next = f_rk4(f, x, u, delT)
%F_RK4 この関数の概要をここに記述
%   詳細説明をここに記述
    k1 = f(x,u);
    k2 = f(x+delT/2 * k1, u);
    k3 = f(x+delT/2 * k2, u);
    k4 = f(x+delT*k3, u);
    x_next = x + delT/6 * (k1 + 2*k2 + 2*k3 + k4);
end