function [vl, ml, xe] = udl(x1, x2, w0, ri, mi)
% This function is for UDL but better

L = x2 - x1;  % span of udl

ra = ri;
%u = 2;

x = linspace(0,L, 10000);   % just making the x coordinates to use

% Plots for udl
    subplot(3, 1, 1) % plot distribution
    hold on
    
    plot([x1, x2], [w0, w0], 'g')
    plot([x1, x1], [0, w0], 'g')
    plot([x2, x2], [0, w0], 'g')

    v = - w0*x + ra;    % force
    m = -w0*x.*x/2 + ra*x + mi; % moment
    
    %d = m.*x.^2/2 - m*L/2.*x - di;
    
    subplot(3, 1, 2)    % plot for shear force
    hold on
    
    plot(linspace(x1, x2,10000), v, 'r')
    
    subplot(3, 1, 3)     % plot for moment
    hold on
    
    %m2 = ra*(x) - w0*x.*x/2 + mr;
    
    plot(linspace(x1, x2, 10000), m, 'b')
    
    % plot for deflection
    %subplot(2, 2, 3)
    %hold on
    %plot(linspace(x1, x2, 10000), d, 'm')
    
    vl = v(10000);
    ml = m(10000);
    xe = x2;
    %de = d(10000);

end
