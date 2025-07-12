%Runge Kutta implementation for solving a system of first degree
%differential equations of the form
%y'i=fi(yi,xi) with Intial Conditions
function [yi] = runge_kutta4sys(f, yi, xi, n, h)
    for j = 1:n
        k1 = h * f(xi(j), yi(j,:)');
        k2 = h * f(xi(j) + h/2, yi(j,:)' + k1/2);
        k3 = h * f(xi(j) + h/2, yi(j,:)' + k2/2);
        k4 = h * f(xi(j) + h, yi(j,:)' + k3);
        yi(j+1,:) = yi(j,:) + (1/6) * (k1' + 2*k2' + 2*k3' + k4');
    end
end