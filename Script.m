%% Coded for the fulfilment of Master's Degree at Politecnico Di Milano
% Author:: Eshwar Bhargav Bhupanam
% Course:: Modeling and Simulation of Aerospace Systems
% Topic:: Non Linear ODE Integration - Runge Kutta, Backward Interpolation
% Year:: 2020-2021

%% Initializing....
clear vars; clc; close all

%% Runge Kutta Method
clear vars; close all; clc;     % Clean variables from the prev. section
f = @(x,t) x-t.^2+1;
x0 = 1/2;
t = [0 2];
CASES = {0.5 'h1';0.2 'h2';0.05 'h3'; 0.01 'h4'};

% Analytical method
tic
[~,x_analytic,h_analytic] = analytic(CASES,t);
CPU_time.analytic = toc;

% Runge Kutta Method
rkorder = 2;                                  % RK2 method
[~,x_rk2,h_rk2,CPU_time_rk2] = rungekutta(CASES,t,x0,f,rkorder);

rkorder = 4;                                  % RK4 method
[tstep,x_rk4,h_rk4,CPU_time_rk4] = rungekutta(CASES,t,x0,f,rkorder);

% Integration error
figure ()
for i = 1:size(CASES,1)
    error_rk2.(CASES{i,2}) = abs(x_analytic.(CASES{i,2})-x_rk2.(CASES{i,2}));
    error_rk4.(CASES{i,2}) = abs(x_analytic.(CASES{i,2})-x_rk4.(CASES{i,2}));
    
    hold on; grid on; legend('show','Location','best');
    plot(tstep.(CASES{i,2}), error_rk2.(CASES{i,2}),'--','DisplayName',strjoin({'RK2 Error @ h =',num2str(CASES{i,1})}),'LineWidth',1.0);
    plot(tstep.(CASES{i,2}), error_rk4.(CASES{i,2}),'-','DisplayName',strjoin({'RK4 Error @ h =',num2str(CASES{i,1})}),'LineWidth',1.5);
end
xlabel ( 'Timespan', 'Interpreter', 'latex' )
ylabel ( 'Integration error \(e(t)\)', 'Interpreter', 'latex')

% Plotting
figure ()
for i = 1:size(CASES,1)
    hold on; grid on; legend('show','Location','best');
    plot(tstep.(CASES{i,2}), x_analytic.(CASES{i,2}),'-','DisplayName',strjoin({'Analytic @ h =',num2str(CASES{i,1})}),'LineWidth',2.0);
    plot(tstep.(CASES{i,2}), x_rk2.(CASES{i,2}),'--','DisplayName',strjoin({'RK2 @ h =',num2str(CASES{i,1})}),'LineWidth',1.0);
end
xlabel ( 'Timespan', 'Interpreter', 'latex' )
ylabel ( 'Solution \(x(t)\)', 'Interpreter', 'latex')

% Plotting
figure ()
for i = 1:size(CASES,1)
    hold on; grid on; legend('show','Location','best');
    plot(tstep.(CASES{i,2}), x_analytic.(CASES{i,2}),'-','DisplayName',strjoin({'Analytic @ h =',num2str(CASES{i,1})}),'LineWidth',2.0);
    plot(tstep.(CASES{i,2}), x_rk4.(CASES{i,2}),'-.','DisplayName',strjoin({'RK4 @ h =', num2str(CASES{i,1})}),'LineWidth',1.0);
end
xlabel ( 'Timespan', 'Interpreter', 'latex' )
ylabel ( 'Solution \(x(t)\)', 'Interpreter', 'latex')

%% Two Dimensional System - Runge Kutta w.r.t Taylor series method ~ H-LAMDA plot
clear vars; close all; clc;     % Clean variables from the prev. section
tol = 1e-2;                     % Tolerance
alpha = 0:pi/180:pi;            % for alpha = [0 pi]
h=[0 5];                        % Initial guess

% Runge Kutta Method
rk_order = 2;                                  % RK2 method
[F.rk2,h_pi.rk2,h_approx.rk2,lamda.rk2] = twodimensionsystem(rk_order,tol, alpha,h);
xlabel ( 'Re\(({h\lambda})\)', 'Interpreter', 'latex' )
ylabel ( 'Imag\(({h\lambda})\)', 'Interpreter', 'latex')

rk_order = 4;                                  % RK4 method
[F.rk4,h_pi.rk4,h_approx.rk4,lamda.rk4] = twodimensionsystem(rk_order,tol, alpha,h);
xlabel ( 'Re\(({h_i\lambda})\)', 'Interpreter', 'latex' )
ylabel ( 'Imag\(({h_i\lambda})\)', 'Interpreter', 'latex')

% Back Interpolation method - BI2 ~ H-ALPHA-THETA plot
theta = {0.4 2 'theta1' 'red';0.1 1.7 'theta2' 'green';0.3 3 'theta3' 'blue';0.7 3 'theta4' 'black';0.9 4 'theta5' 'cyan'};
figure()
labels = {};
for k = 1 : size(theta,1)
    [B,h_alpha.(theta{k,3}),lamda.(theta{k,3})] = plotBI2stability(theta{k,1},alpha,theta{k,2},theta{k, 4});
    labels{end+1} = strjoin({strcat('\color{',theta{k, 4},'}'),strjoin({'\theta =',num2str(theta{k,1})})});
end
legend(labels)
xlabel ( 'Re\(({h\lambda})\)', 'Interpreter', 'latex' )
ylabel ( 'Imag\(({h\lambda})\)', 'Interpreter', 'latex')
Operator = B;

%% Runge Kutta and Backward Interpolation Method
clearvars; close all; clc;     % Clean variables from the prev. section

% Initial Value Problem using RK4
h = 0.1; % h = (b-a)/N (Input)
B = [-180.5, 219.5; 179.5, -220.5]; % Input
IVC = [1;1]; %Initial condition

a = 0; % Lower bound
b = 5; % Upper bound
N = (b-a)/h;
teval=0:h:5;

% From ---> https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
% In relation with the Taylor series and RK4
% for xdot = B*x
rk_order = 4;
[x_rk4,~] = RK_Taylor(rk_order,IVC,B,h,N);

% Initial Value Problem using BI2_0.1
theta = 0.1;
[x_BI2,~] = BI2(theta,IVC,B,N,h);

% Analytical solution using feval
x_analytic = fun_analytic(B,N,teval,IVC);

error_rk4 = abs(x_analytic-x_rk4);
error_BI2 = abs(x_analytic-x_BI2);

figure(1);hold on; grid on;
plot(teval,x_analytic,'b-x', teval,x_rk4,'k-o');
xlabel ( 'Time(t) \([s]\)', 'Interpreter', 'latex' )
ylabel ( 'Solution $x(t)$', 'Interpreter', 'latex')
hold off;

figure(2);hold on; grid on;
plot(teval,x_analytic,'b-x', teval,x_BI2,'r-^');
xlabel ( 'Time(t) \([s]\)', 'Interpreter', 'latex' )
ylabel ( 'Solution $x(t)$', 'Interpreter', 'latex')
hold off;

figure(3);hold on; grid on;
plot(teval,error_rk4,'k-x',teval,error_BI2,'r-^');
hold off;

% Numerical Stability plot
lamda = eig(B);
% https://www.mathworks.com/matlabcentral/answers/366363-error-on-contour-
% Plot-stability-region-runngekutta
figure(4); hold on; grid on;
[X,Y] = meshgrid(-5:0.01:5,-5:0.01:5);
Mu = X+1i*Y;
R = 1 + Mu + .5*Mu.^2 + (1/6)*Mu.^3 + (1/24)*Mu.^4;
Rhat = abs(R);
contour(X,Y,Rhat,[1 1]);
plot(real(h*lamda),imag(h*lamda),'k.','Markersize',30,'DisplayName','HI');
xlabel ( 'Re\(({h\lambda})\)', 'Interpreter', 'latex' )
ylabel ( 'Imag\(({h\lambda})\)', 'Interpreter', 'latex')

% Plot-stability-region-BI2
figure(5)
theta = {0.1 1.7, 'red'};
alpha = 0:pi/180:pi;
labels={};
for k = 1:size(theta,1)
    [B,h_alpha,eigA] = plotBI2stability(theta{k,1},alpha,theta{k,2}, theta{k,3});
    labels{end+1} = strjoin({strcat('\color{',theta{k, 3},'}'),strjoin({'\theta =',num2str(theta{k,1})})});
end
legend(labels)
h = 0.1; hold on;
plot(real(h*lamda),imag(h*lamda),'k.','Markersize',30,'DisplayName',strjoin({'h @',num2str(h)}));
xlabel ( 'Re\(({h\lambda})\)', 'Interpreter', 'latex' )
ylabel ( 'Imag\(({h\lambda})\)', 'Interpreter', 'latex')

%% Functions

function [tstep,x,h] = analytic(CASES,t)
for i = 1:size(CASES,1)
    h = cell2mat(CASES(i,1));
    tstep.(CASES{i,2}) = t(1):h:t(end);
    x.(CASES{i,2}) = zeros(max(length(tstep.(CASES{i,2}))),1);
    x.(CASES{i,2}) = tstep.(CASES{i,2}).^2 +2.*tstep.(CASES{i,2}) +1 -(1/2).*exp(tstep.(CASES{i,2}));
end
end

function [tstep,x,h,CPU_time] = rungekutta(CASES,t,x0,f,rkorder)

for i = 1:size(CASES)
    h = cell2mat(CASES(i,1));
    tstep.(CASES{i,2}) = t(1):h:t(end);
    x.(CASES{i,2}) = zeros(1,max(length(tstep.(CASES{i,2}))));
    xpre.(CASES{i,2}) = zeros(max(length(tstep.(CASES{i,2})))-1,1);
    x.(CASES{i,2})(1,1)=x0;
    for j = 1:length(tstep.(CASES{i,2}))-1
        if rkorder == 2
            tic
            alpha = [1;1];
            beta = [1 0;1/2 1/2];
            % Predictor
            xpre.(CASES{i,2})(j) = x.(CASES{i,2})(j) + ...
                beta(1,1)*h*feval(f,x.(CASES{i,2})(j),tstep.(CASES{i,2})(j));
            % Corrector
            x.(CASES{i,2})(j+1) = x.(CASES{i,2})(j) + alpha(2)*h*(beta(2,1)*...
                (feval(f,x.(CASES{i,2})(j),tstep.(CASES{i,2})(j)))+(beta(2,2)*...
                (feval(f,xpre.(CASES{i,2})(j),(tstep.(CASES{i,2})(j)+(alpha(1)*h))))));
            CPU_time.(CASES{i,2}) = toc;
        elseif rkorder == 4
            tic
            alpha = [1/2;1/2;1;1];
            beta = [1/2 0 0 0;...
                0 1/2 0 0;...
                0 0 1 0;...
                1/6 1/3 1/3 1/6];
            % Predictor
            xpre1.(CASES{i,2})(j) = x.(CASES{i,2})(j) + ...
                beta(1,1)*h*feval(f,x.(CASES{i,2})(j),tstep.(CASES{i,2})(j));
            xpre2.(CASES{i,2})(j) = x.(CASES{i,2})(j) + ...
                beta(2,2)*h*feval(f,xpre1.(CASES{i,2})(j),tstep.(CASES{i,2})(j)+(alpha(1)*h));
            xpre3.(CASES{i,2})(j) = x.(CASES{i,2})(j) + ...
                beta(3,3)*h*feval(f,xpre2.(CASES{i,2})(j),tstep.(CASES{i,2})(j)+(alpha(2)*h));
            % Corrector
            x.(CASES{i,2})(j+1) = x.(CASES{i,2})(j) +...
                alpha(4)*h*(beta(4,1)*(feval(f,x.(CASES{i,2})(j),tstep.(CASES{i,2})(j)))+...
                (beta(4,2)*(feval(f,xpre1.(CASES{i,2})(j),(tstep.(CASES{i,2})(j)+(alpha(1)*h)))))+...
                (beta(4,3)*(feval(f,xpre2.(CASES{i,2})(j),(tstep.(CASES{i,2})(j)+(alpha(2)*h)))))+...
                (beta(4,4)*(feval(f,xpre3.(CASES{i,2})(j),(tstep.(CASES{i,2})(j)+(alpha(3)*h))))));
            CPU_time.(CASES{i,2}) = toc;
        end
        
    end
end
end

function [x_rk4,term_rk4] = RK_Taylor(rk_order,IVC,B,h,N)
x_rk4 = zeros(length(B), N+1);
x_rk4(:,1) = IVC;
expansion_terms = 1;
for k=1:rk_order
    expansion_terms = (expansion_terms)*((B.*h)/k);
    if k == 1
        term_rk4 = eye(length(B))+ expansion_terms;
    else
        term_rk4 = term_rk4 + expansion_terms;
    end
    for i = 1:N
        x_rk4(:,i+1)= term_rk4*x_rk4(:,i);
    end
end
end

function [h_approx,eigA] = find_hRK(alpha,h,rk_order,tol)
h_approx = zeros(1,length(alpha));
eigA = zeros(2,length(alpha));
figure()
for i=1:length(alpha)
    A = find_A(alpha(i));
    eigA(:,i) = eig(A);
    condition = 1+tol;
    iter = 0;
    hnew = h;
    while abs(condition) > tol
        iter = iter+1;
        h_guess = (h(1)+h(2))/2;
        [~,F] = RK_Taylor(rk_order,0,A,h_guess,0);
        condition = -1+max(abs(eig(F)));
        if condition > 0
            h(2) = h_guess;
        else
            h(1) = h_guess;
        end
    end
    h_approx(1,i) = h_guess;
    clear h; h = hnew;
    hold on; grid on;
    % Plotting h-lamda plane
    plot(real(eigA(:,i)*h_approx(1,i)),imag(eigA(:,i)*h_approx(1,i)),'r*')
    legend ({'h @ \alpha = [0 \pi]'})
    if rk_order == 4
        hold on;
        CASES = {0.5 'h1';0.2 'h2';0.05 'h3'; 0.01 'h4'};
        for j = 1:size(CASES,1)
            hi = CASES{j,1};
            plot(real(hi),imag(hi),'.','MarkerSize',25,'DisplayName', strjoin({'h @',num2str(CASES{j,1})}))
        end
    end
end
hold on;
plot(real(eigA(:,end)*h_approx(1,end)),imag(eigA(:,end)*h_approx(1,end)),...
    'g.','MarkerSize',20,'DisplayName','h @ \alpha = \pi')
end

function [F,h_pi,h_approx,lamda] = twodimensionsystem(rk_order,tol,alpha,h)
[h_approx,lamda] = find_hRK(alpha,h,rk_order,tol);
h_pi = h_approx(end);      % for alpha = pi "hrk2_pi" is the solution
syms h A
[~,F] = RK_Taylor(rk_order,0,A,h,0);
end

function [x_BI2,term_BI2] = BI2(theta,IVC,B,N,h)
x_BI2(:,1) = IVC;
P = eye(length(B))-((1-theta)*h*B);
Q = eye(length(B))+(theta*h*B);
term_BI2 = mldivide(P,Q);
for i = 1:N
    x_BI2(:,i+1) = term_BI2*x_BI2(:,i);
end
end

function [B,h_alpha,eigA] = plotBI2stability(theta,alpha,x,color)
h_alpha = zeros(1,length(alpha));
eigA = zeros(2,length(alpha));
for i = 1:length(alpha)
    A = find_A(alpha(i));
    eigA(:,i) = eig(A);syms h;
    h_temp =@(h) find_hBI(A,h,theta);
    h_alpha(1,i)=fsolve(h_temp,x);
    hold on; grid on;
    plot(real(eigA(:,i)*h_alpha(1,i)),imag(eigA(:,i)*h_alpha(1,i)),strcat(color,'.'));
end
% labels{end+1} = strjoin({'\theta =',num2str(theta)});
% legend(labels)
syms A theta
[~,B] = BI2(theta,0,A,0,h);
end

function [Sol_analytic]=fun_analytic(B,N,t,x0)
% x = exp(B.*t)*x0;
Sol_analytic = zeros(length(x0),N+1);
for j = 1:N+1
    f =@(t) expm(B.*t)*x0;
    Sol_analytic(:,j) = feval(f,t(j));
end
end

function A = find_A(alpha)
A = [0,1; -1,2*cos(alpha)];
end

function h_temp = find_hBI(B,h,theta)
P = eye(length(B))-((1-theta)*h*B);
Q = eye(length(B))+(theta*h*B);
term_BI2 = mldivide(P,Q);

h_temp = -1+max(abs(eig(term_BI2)));
end
