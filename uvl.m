function [ vl, ml, xe ] = uvl(x1, f1, x2, f2, ri, mi )

% This function is for not udl

L = x2 - x1;  % span of uvl
x = linspace(0,L, 10000);
w0 = abs(f2 - f1); % difference in loads

if f2>f1
    k = 1;
else
    k = -1;
end

c = ( f1 + k*w0/L*x + 2*f1)./(3*(f1 + f1 + k*w0/L*x));

v = -1/2 * (f1 + f1 + k*w0*x/L).*x + ri; % shear force
m =  -1/2*(f1 + f1 + k*w0/L*x).*x.*c.*x + ri*x + mi ; % moment

subplot(3, 1, 1) % plot distribution of forces
hold on

plot([x1, x1], [0, f1], 'black')
plot([x1, x2], [f1, f2], 'black')
plot([x2, x2], [0, f2], 'black')

subplot(3, 1, 2) % plot shear force
hold on

plot(linspace(x1, x2, 10000), v, 'r')

subplot(3, 1, 3) % plot moment
hold on

plot(linspace(x1, x2, 10000), m, 'b')

xe = x2;
vl = v(10000);
ml = m(10000);


end

