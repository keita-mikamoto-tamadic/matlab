function x_dec = x2rad(theta, phi, a, l)
%X2RAD この関数の概要をここに記述
    x_wheel = a * (theta + phi);
    x_pendulam = l * sin(theta) + x_wheel;
    y_pendulam = l * cos(theta);
    x_phi = a * sin(theta + phi) +x_wheel;
    y_phi = a * cos(theta + phi);
    x_dec = [x_wheel; x_pendulam; y_pendulam; x_phi; y_phi];
end

