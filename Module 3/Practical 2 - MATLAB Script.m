clc, clearvars;

dydt = @(t, y) (y + 1);

t_start = 0.01;
t_final = 1;
y_final = 0.1;
h = 0.01; 

tspan = t_final : -h : t_start;

% 1(a)

y_euler = zeros(1, length(tspan));
y_euler(1) = y_final;
for i = 1:length(tspan)-1
    y_euler(i+1) = y_euler(i) - h * dydt(i, y_euler(i));
end

[t_runge_kutta, y_runge_kutta] = ode45(dydt, tspan, y_final);

tspan = fliplr(tspan);  
y_euler = fliplr(y_euler);
y_runge_kutta = flipud(y_runge_kutta);

figure;
plot(tspan, y_euler, 'b', 'DisplayName', 'Euler Method');
hold on;
plot(tspan, y_runge_kutta, 'r', 'DisplayName', 'Runge-Kutta Method');
xlabel('Time (s)');
ylabel('y(t)');
title('Comparison of Euler and Runge-Kutta Methods');
legend('Location', 'Best');
grid on;


% 1(b)

h_values = [0.001, 0.01, 0.05, 0.1, 0.15, 0.2]; 

figure;

for i = 1:length(h_values)

    h = h_values(i);

    tspan = t_final : -h : t_start;

    if tspan(end) ~= t_start
        tspan = [tspan, t_start];  
    end
    
    y_euler = zeros(1, length(tspan));
    y_euler(1) = y_final;
    for j = 1:length(tspan)-1
        y_euler(j+1) = y_euler(j) - h * dydt(j, y_euler(j));
    end

    [t_runge_kutta, y_runge_kutta] = ode45(dydt, tspan, y_final);
    
    tspan = fliplr(tspan);
    y_euler = fliplr(y_euler);
    y_runge_kutta = flipud(y_runge_kutta);
    
    subplot(3, 2, i);
    plot(tspan, y_euler, 'b', 'DisplayName', 'Euler Method');
    hold on;
    plot(tspan, y_runge_kutta, 'r', 'DisplayName', 'Runge-Kutta Method');
    xlabel('Time (s)');
    ylabel('y(t)');
    title(['h = ', num2str(h)]);
    legend('Location', 'Best');
    grid on;

end

% 2(a)

tspan=0:0.01:100;
h = 0.01;

x_euler=zeros(13,length(tspan));
x_euler(1,1) = 0.5;
x_euler(2,1) = 3.0E-4;
x_euler(3,1) = 1.0E-4;
x_euler(4,1) = 3.4;
x_euler(5,1) = 0.06;
x_euler(6,1) = 0.96;
x_euler(7,1) = 0.03;
x_euler(8,1) = 2.464;
x_euler(9,1) = 0.04;
x_euler(10,1) = 0.08528;
x_euler(11,1) = 0.408;
x_euler(12,1) = 0.114;
x_euler(13,1) = 0.1;

for i = 1:length(tspan)-1
    x_euler(:,i+1) = x_euler(:,i) + h * tca_xdot(i,x_euler(:,i));
end

figure;
plot(tspan,x_euler);
xlabel('Time (min)');
ylabel('Concentration (mM)');
title('TCA Model using Euler with time step = 0.01');
legend('aca','oaa','coa','cit','icit','akg','ssa','suc','sca','fa','mal','gly','biosyn','Location', 'Best');
grid on;


tspan=0:0.001:100;
h = 0.001;

x_euler=zeros(13,length(tspan));
x_euler(1,1) = 0.5;
x_euler(2,1) = 3.0E-4;
x_euler(3,1) = 1.0E-4;
x_euler(4,1) = 3.4;
x_euler(5,1) = 0.06;
x_euler(6,1) = 0.96;
x_euler(7,1) = 0.03;
x_euler(8,1) = 2.464;
x_euler(9,1) = 0.04;
x_euler(10,1) = 0.08528;
x_euler(11,1) = 0.408;
x_euler(12,1) = 0.114;
x_euler(13,1) = 0.1;

for i = 1:length(tspan)-1
    x_euler(:,i+1) = x_euler(:,i) + h * tca_xdot(i,x_euler(:,i));
end

figure;
plot(tspan,x_euler);
xlabel('Time (min)');
ylabel('Concentration (mM)');
title('TCA Model using Euler with time step = 0.001');
legend('aca','oaa','coa','cit','icit','akg','ssa','suc','sca','fa','mal','gly','biosyn','Location', 'Best');
grid on;


figure;
xlabel('Time (min)');
ylabel('Concentration (mM)');
BIOMD0000000219;
title('TCA Model using Biomodels with ode23');
legend('aca','oaa','coa','cit','icit','akg','ssa','suc','sca','fa','mal','gly','biosyn','Location', 'Best');
grid on;
