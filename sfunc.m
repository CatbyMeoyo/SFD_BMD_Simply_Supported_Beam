function [] = sfunc( l, a, b )

if nargin == 3 && strcmp(a, 'udl') && strcmp(b, 'uvl') == true
    u = 2;
    v = 3;
elseif nargin == 2 && strcmp(a, 'uvl') == true
    u = 0;
    v = 3;
elseif nargin == 2 && strcmp(a, 'udl') == true
    u = 2;
    v = 0;
elseif nargin == 1
    u = 0;
    v = u;
end

p = 1;
    
% initialize
ip = 0;
iu = 0;
iv = 0;
p_forces = zeros(l*3, 5);   % originally of shape n by 3
% end initialize

fprintf('Position measured from left end\n')

% take inputs for point loads
fprintf('\n Input Point Loads\n')

    while (true)
        input_forces = input('[position, Applied force, moment] = ');
        if isempty(input_forces) == true
            break % EXTERMNATE!
        else
            ip = ip + 1;
            p_forces(ip, :) = [input_forces, 0, p];
            % making an input matrix
        end
        
        if input_forces(1) > l
            error('Position greater than length!')
        end
    end

p_forces = p_forces(1:ip, :); % cutting the matrix to exclude initialized zeros
% end of taking inputs

% sum forces for point loads
    % total force for point forces
    p_tf = 0;
    % moment from the right end for point forces
    ptf_by_l = 0;
        for  jp = 1 : length(p_forces(:,2))
            p_tf = p_tf + p_forces(jp,2);

            ptf_by_l = ptf_by_l + (l - p_forces(jp,1))*p_forces(jp,2) - p_forces(jp,3);
        end
    
% end of summing up forces for point loads

% take inputs for ud loads
    if u == 2
    fprintf('\n Input uniformly distributed Loads\n')
        u_forces = zeros(l*3, 5);   % originally of shape n by 3
        u_forces_d = zeros(l*3, 5);   % false positions, to correct the plot, originally of shape n by 3
        while (true)
            input_forces2 = input('[start position, end position, force per length] = ');
            if isempty(input_forces2) == true
                break % EXTERMNATE!
            else
                iu = iu + 1;
                u_forces(iu, :) = [input_forces2*(1+ 10^-(3+8)), 0, u]; % multiplying with 1 + 10^(3+8) to correct the plot
                u_forces_d(iu, :) = [input_forces2(2)*(1 + 10^-38), 0, 0, 0, p];
                % making an input matrix
            end
            
            if input_forces2(1) > l 
                error('Position greater than length!')
            elseif input_forces2(2) > l
                error('Position greater than length!')
            end
            
        end

    u_forces = u_forces(1:iu, :); % cutting the matrix to exclude initialized zeros
    u_forces_d = u_forces_d(1:iu, :); % cutting the matrix to exclude initialized zeros
    end
% end of taking inputs

% sum forces for ud loads
if u == 2
    % total force for point forces
    u_tf = 0;
    % moment from the right end for point forces
    utf_by_l = 0;
        for  ju = 1 : length(u_forces(:,3))
            u_tf = u_tf + u_forces(ju,3)*(u_forces(ju,2) - u_forces(ju,1));
            
            utf_by_l = utf_by_l + u_forces(ju,3)*(u_forces(ju,2) - u_forces(ju,1))*(l- u_forces(ju,2) +(u_forces(ju,2) - u_forces(ju,1))/2);
        end
end
    
% end of summing up forces for ud loads


% take inputs for uv loads
    if v == 3
    fprintf('\n Input uniformly varying Loads\n')
        v_forces = zeros(l*3, 5);   % originally of shape n by 4
        v_forces_d = zeros(l*3, 5);   % false positions, to correct the plot, originally of shape n by 4
        while (true)
            input_forces3 = input('[start position, start unit force, end position, end unit force] = ');
            if isempty(input_forces3) == true
                break % EXTERMNATE!
            else
                iv = iv + 1;
                v_forces(iv, :) = [input_forces3*(1+10^-(3+8)), v]; % multiplying with 1 + 10^(3+8) to correct the plot
                v_forces_d(iv, :) = [input_forces3(3)*(1+10^-38), 0, 0, 0, p];
                % making an input matrix
            end
            
            if input_forces3(1) > l 
                error('Position greater than length!')
            elseif input_forces3(3) > l
                error('Position greater than length!')
            end
            
        end

    v_forces = v_forces(1:iv, :); % cutting the matrix to exclude initialized zeros
    v_forces_d = v_forces_d(1:iv, :); % cutting the matrix to exclude initialized zeros
    end
% end of taking inputs

%  sum forces for uv loads

if v == 3
    % total force for point forces
        v_tf = 0;
    % moment from the right end for point forces
        vtf_by_l = 0;
        
        for  jv = 1 : length(v_forces(:,1))
            v_tf = v_tf + 1/2 * (v_forces(jv, 2) + v_forces(jv,4)) * (v_forces(jv,3) - v_forces(jv, 1));
            
            vtf_by_l = vtf_by_l +  (1/2 * (v_forces(jv, 2) + v_forces(jv,4)) * (v_forces(jv,3) - v_forces(jv, 1))) * ((v_forces(jv,3) - v_forces(jv, 1))*(v_forces(jv,4) + 2*v_forces(jv,2))./(3*(v_forces(jv,2) + v_forces(jv,4))) +l - v_forces(jv,3));
        end
end

% end of summing up forces for ud loads

% now to arange and sort, somewhat

% note:
% tdf = total acting loads
% ra = reaction at left most end
% rb = reaction at right most end

zero_pointl = [0, 0, 0, 0, 1]; % at point zero, the forces

if u==2 && v==3
    tforces_1 = [zero_pointl; p_forces; u_forces; v_forces];
    
    tdf = p_tf + u_tf + v_tf;
    ra = (ptf_by_l + utf_by_l + vtf_by_l)/l;
    rb = tdf - ra;
    
    tforces_2 = [ tforces_1; u_forces_d ;v_forces_d];
    
elseif u == 2 && v == 0
    tforces_1 = [zero_pointl; p_forces; u_forces];
    
    tdf = p_tf + u_tf;
    ra = (ptf_by_l + utf_by_l)/l;
    rb = tdf - ra;
    
    tforces_2 = [ tforces_1; u_forces_d];
    
elseif u == 0 && v == 3
    tforces_1 = [zero_pointl; p_forces; v_forces];
    
    tdf = p_tf + v_tf;
    ra = (ptf_by_l + vtf_by_l)/l;
    rb = tdf - ra;
    
    tforces_2 = [ tforces_1; v_forces_d];
    
else
    tforces_1 = [zero_pointl; p_forces];
    
    tdf = p_tf;
    ra = ptf_by_l/l;
    rb = tdf - ra;
    
    tforces_2 = tforces_1;
    
end

force_mat = sortrows(tforces_2, [1, 5]); % matrix containing all the forces

% now, let's plot...

i = length(force_mat(:, 1)); % total iteration number

subplot(3, 1, 1) % plots of loads are here
    hold on
    
    xlim([0, 1.1*l])
    plot([0,l], [0,0], 'm')
    
    xlabel('Distance')
    ylabel('Magnitude of Loads')
    title('Load Distribution')

subplot(3, 1, 2) % plots for shear forces are here
    
    hold on
    xlim([0, 1.1*l])
    %ylim([-1.5*max(abs(ra)), 1.5*max(abs(ra))])
    plot([0,l], [0,0], 'b')
    
    xlabel('Distance')
    ylabel('Shear Force')
    title('Shear Force Diagram')

       
subplot(3, 1, 3) % plots for moments are here
    
    hold on
    xlim([0, 1.1*l])
    plot([0,l], [0,0], 'b')
    
    xlabel('Distance')
    ylabel('Moment')
    title('Bendng Moment Diagram')

ri = ra; % initial reaction force for every point
mi = 0; % initial moment for every point
xp = 0;

xps = [force_mat(:, 1); l];

for ti = 1:i 
    if force_mat(ti, 5) == 1 % point loads
        iforce = force_mat(ti, 1:3);
        xl = xps(ti+1);
        [ri, mi, xp] = pointl(xl, xp, ri, mi, iforce); % plotting and also updating the values
    end
    
    if force_mat(ti, 5) == 2
        xu1 = force_mat(ti, 1);
        xu2 = force_mat(ti, 2);
        w0 = force_mat(ti, 3);
        [ri, mi, xp] = udl(xu1, xu2, w0, ri, mi); % plotting and also updating the values
    end
    
    if force_mat(ti, 5) == 3
        xv1 = force_mat(ti, 1);
        xv2 = force_mat(ti, 3);
        fv1 = force_mat(ti, 2);
        fv2 = force_mat(ti, 4);
        [ri, mi, xp] = uvl(xv1, fv1, xv2, fv2, ri, mi);
    end
end

% residue plot
    subplot(3, 1, 1)
    hold off
    
    subplot(3, 1, 2)
    
    plot([xp,l], [ri, -rb], 'r')
    plot([l,l], [0, -rb], 'r')
    hold off
    
    subplot(3, 1, 3)
    
    plot([xp,l], [mi, 0], 'b')
    hold off

end

