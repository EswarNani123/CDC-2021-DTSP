function fin_head = final_heading(p1,p2,rho)
ini_head = mod(p1(3),2*pi);
x = (p2(1)-p1(1))*cos(-ini_head)-(p2(2)-p1(2))*sin(-ini_head);
y = (p2(2)-p1(2))*cos(-ini_head)+(p2(1)-p1(1))*sin(-ini_head);
dc_q = sqrt(x^2+(abs(y)-rho)^2);
theta_c = atan2(x,rho-abs(y));
df_q = sqrt(x^2+(abs(y)+rho)^2);
alpha = acos((5*rho^2-df_q^2)/(4*rho^2));
if sqrt(x^2+(y-rho)^2) < rho % RL
    fin_head = -asin((dc_q*sin(theta_c))/df_q)-asin((rho*sin(alpha))/df_q)+2*pi-alpha+ini_head;
end
if sqrt(x^2+(y+rho)^2) < rho % LR
    fin_head = asin((dc_q*sin(theta_c))/df_q)+asin((rho*sin(alpha))/df_q)-2*pi+alpha+ini_head;
end
if (sqrt(x^2+(y-rho)^2) > rho) && (sqrt(x^2+(y+rho)^2) > rho) && (y > 0)
    fin_head = theta_c-acos(rho/dc_q)+ini_head; % LS
end
if ((sqrt(x^2+(y-rho)^2) > rho) && (sqrt(x^2+(y+rho)^2) > rho) && (y < 0)) || ((sqrt(x^2+(y-rho)^2) > rho) && (sqrt(x^2+(y+rho)^2) > rho) && (y == 0) && x < 0)
    fin_head = -theta_c+acos(rho/dc_q)+ini_head; % RS
end
if sqrt(x^2+(y-rho)^2) == rho 
    fin_head = theta_c-acos(rho/dc_q)+ini_head; % L
end
if sqrt(x^2+(y+rho)^2) == rho
    fin_head = -theta_c+acos(rho/dc_q)+ini_head; % R
end
if (sqrt(x^2+(y-rho)^2) > rho) && (sqrt(x^2+(y+rho)^2) > rho) && (y == 0) && (x > 0)
    fin_head = ini_head; % S
end
