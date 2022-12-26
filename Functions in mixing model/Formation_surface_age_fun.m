function [F,S,m_tp,m,Krw_real,s] = Formation_surface_age_fun(t,Mud,Mdd,Mc,Krw)
% Formation_surface_age_fun.m calculate formation age distribution
% and surface age distribution simultaneously
%
% Meng Guo, Yale University
% Summer, 2019

t = t';
Mcp = Mc(end);% present-day continental crust mass
dt = t(2)-t(1);% calcualte timestep dt

nt=length(t);
m = nan(nt,nt); % m(tau,t): formation age distribution at time t
s = nan(nt,nt); % s(tau,t): surface age distribution at time t

F = nan(size(t));
S = nan(size(t));
m(:,1) = zeros(1,nt); % when t=0
s(:,1) = zeros(1,nt);


delta = zeros(size(t)); 
delta = delta';       

Krw_real= zeros(nt,1);

for i = 2:nt 
  delta(:) = 0;
  delta(i) = 1/dt; % because time is discretized with dt
  % note: dt*delta = 1 in the following
  
  if Mc(i)==0
    m(:,i) = m(:,i-1) + dt*Mud(i)*delta; 
    s(:,i) = s(:,i-1) + dt*Mud(i)*delta;
    % when there is no continental crust one timestep before, 
    % there is only generation of new crust, no recycling
  else
    % formation age distribution
    m(:,i) = m(:,i-1) + dt*(Mud(i)*delta - (Mdd(i)/Mc(i))*m(:,i-1));
    for j = 1:nt
        if (m(j,i)<0)
           m(j,i) = 0;
        end 
    end 

    % surface age distribution
    % first consider the effect of recycling
    s(:,i) = s(:,i-1) + dt*(Mud(i)*delta - (Mdd(i)/Mc(i))*s(:,i-1));
    
    % then consider the effect of reworking
    tmpval(i) = Krw(i)/Mc(i)*dt;
    if tmpval(i)  > 1
        tmpval(i) = 1;
        s(1:(i-1),i) = 0;
    else
        tmpval(i) = Krw(i)/Mc(i)*dt;
        s(1:(i-1),i) = s(1:(i-1),i)*(1-tmpval(i));
    end
    sum_val = sum(s(1:(i-1),i)*tmpval(i));
    s(i,i) = s(i,i)+sum_val;
    Krw_real(i) = tmpval(i)*Mc(i)/dt;
  end
  for j = 1:nt;
    if (s(j,i)<0)
      s(j,i) = 0;
    end
  end

  if 0 % set this to 1 to see intermedieate steps
    if mod(i,100) == 0
      figure(10);
      %plot(t,m(:,i),'r-',t,s(:,i),'b--');
      M = cumsum(m(1:i,i))*dt; M = M/Mcp;
      S = cumsum(s(1:i,i))*dt; S = S/Mcp;
      plot(t(1:i),M,'b-',t(1:i),S,'r--'); 
      title(['t=' num2str(t(i))]);
      %    axis([0 max(t) 0 1.01]);
      %    pause;
      pause(0.01);
    end
  end
end

m_tp = m(:,nt);% m(tau,tp) is the present-day formation age distribution 
s_tp = s(:,nt);
F(1) = m_tp(1);
S(1) = s_tp(1);

% calculate the cumulative formation age and surface age distributions
for j = 2:nt % in tau domain
  F(j) = F(j-1) + m_tp(j)*dt; 
  S(j) = S(j-1) + s_tp(j)*dt;
end

% normalize the cumulative formation age distribution according to present-day
% continental crust mass
F = F/Mcp;  
F = F';

% normalize the cumulative surface age distribution according to present-day
% continental crust mass
S = S/Mcp;
S = S';




