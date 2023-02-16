function [ vl, ml, xe ] = pointl(xl, xp, r, mr, iforce)

% This function is for point loads

forces = [[xp, -r, mr]; iforce];

       % plot the forces
       subplot(3, 1, 1)
       hold on
       
       plot([iforce(1), iforce(1)], [0, iforce(2)], 'r') % plot shear force
       plot([iforce(1), iforce(1)], [0, iforce(3)], 'b') % plot moment

% some initialization is done here
i = length(forces(:,1));
x = [forces(:,1); xl]; % matrix containing positions
f = forces(:,2); % matrix containing forces
m0s = zeros(1, i+2);
m0 = zeros(1, i+2);
m = [forces(:,3); 0];
mcn = zeros(1, 10000);
%d = zeros(1, 10000);

% plot-a-lot!
fcn = 0;

   for j = 1:length(x) - 1
       xn = x(j);
       xm = x(j + 1);
       fcn = fcn - f(j); % shear force
       m0s(1, j) = mcn(1, 10000);
       mcn = fcn*linspace(xn, xm, 10000) + sum(f(1:j).*x(1:j)) + sum(m(1:j)); % bending moment
       m0(1, j) = mcn(1, 1);
       mn = m0s(j);
       mm = m0(j);
       
       %xd = linspace(xn, xm, 10000);
       %d = mcn.*xd.^2/2 + di - trapz(xd, mcn)/(xm - mn).*xd;
              
       % plot shear force diagram
       subplot(3, 1, 2)
       hold on
       
       plot([xn, xm], [fcn, fcn], 'r') % plot the horizontal lines
       plot([xn,xn], [fcn, fcn + f(j)],'r') % plot the vertical lines
       
       % plot bending moment diagram
       subplot(3, 1, 3)
       hold on
       
       plot(linspace(xn, xm, 10000), mcn, 'b')
       plot([xn, xn], [mn, mm], 'b')
       
       %plot([l,l], [0, me], 'g')
       
       %subplot(2, 2, 3)
       %hold on
       %plot(linspace(xn, xm, 10000), -d, 'b')
       
   end
vl = fcn;
ml = mcn(10000);
xe = xm;
%de = d(10000);
   
end

