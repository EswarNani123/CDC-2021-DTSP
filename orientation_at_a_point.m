function alpha = orientation_at_a_point(a,b,c) % function takes three consecutive points a,b and c from 
% ETSP tour sequence
theta_a = atan2(a(2)-b(2),a(1)-b(1)); % Angle of a w.r.t b
theta_c = atan2(c(2)-b(2),c(1)-b(1)); % Angle of c w.r.t b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We do not want negative angles, so the following code convert it to positive 
if theta_a < 0
    theta_a = theta_a+2*pi;
end
if theta_c < 0
    theta_c = theta_c+2*pi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following code rotates the point c
if theta_a > theta_c
    x = c(1)-b(1); % Subtracting point b from c
    y = c(2)-b(2);
    P = [x;y];
    theta = abs(theta_a-theta_c)/2;
    rot_matrix = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    Q = rot_matrix*P;
    Q = Q+b'; % Adding back point b
    Q_minus = [0 1;-1 0]*[Q(1)-b(1);Q(2)-b(2)];
    Q_minus = Q_minus+b';
    theta_Q_minus = atan2(Q_minus(2)-b(2),Q_minus(1)-b(1)); % Angle of Q_minus w.r.t b
    if theta_Q_minus < 0 % Making sure angle is positive
        theta_Q_minus = theta_Q_minus+2*pi;
    end
    alpha = theta_Q_minus;
else
    x = c(1,1)-b(1,1); % Subtracting point b from c
    y = c(1,2)-b(1,2);
    P = [x;y];
    theta = -abs(theta_a-theta_c)/2;
    rot_matrix = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    Q = rot_matrix*P;
    Q = Q+b'; % Adding back point b
    Q_plus = [0 -1;1 0]*[Q(1)-b(1);Q(2)-b(2)]; % Rotates Q by pi/2 
    Q_plus = Q_plus+b';
    theta_Q_plus = atan2(Q_plus(2)-b(2),Q_plus(1)-b(1)); % Angle of Q_plus w.r.t b
    if theta_Q_plus < 0 % Making sure angle is positive
    theta_Q_plus = theta_Q_plus+2*pi;
    end
    alpha = theta_Q_plus;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Q_plus = [0 -1;1 0]*[Q(1)-b(1);Q(2)-b(2)]; % Rotates Q by pi/2 
% Q_plus = Q_plus+b';
% theta_Q_plus = atan2(Q_plus(2)-b(2),Q_plus(1)-b(1)); % Angle of Q_plus w.r.t b
% if theta_Q_plus < 0 % Making sure angle is positive
%     theta_Q_plus = theta_Q_plus+2*pi;
% end
% Q_minus = [0 1;-1 0]*[Q(1)-b(1);Q(2)-b(2)];
% Q_minus = Q_minus+b';
% theta_Q_minus = atan2(Q_minus(2)-b(2),Q_minus(1)-b(1)); % Angle of Q_minus w.r.t b
% if theta_Q_minus < 0 % Making sure angle is positive
%     theta_Q_minus = theta_Q_minus+2*pi;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Choosing the angle of Q_plus or Q_mins which ever is nearest to c
% if abs(theta_Q_plus-theta_c)>abs(theta_Q_minus-theta_c) || abs(theta_Q_plus-theta_c)>abs(2*pi-(theta_Q_minus-theta_c))
%     alpha = theta_Q_minus;
% else
%     alpha = theta_Q_plus;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end