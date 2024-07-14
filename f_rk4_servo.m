function x_next = f_rk4_servo(f, x, fd, ref, delT)
%F_RK4 この関数の概要をここに記述
%   詳細説明をここに記述
    k1 = f(x, fd, ref);
    k2 = f(x+delT/2 * k1, fd, ref);
    k3 = f(x+delT/2 * k2, fd, ref);
    k4 = f(x+delT * k3, fd, ref);
    x_next = x + delT/6 * (k1 + 2*k2 + 2*k3 + k4);
end