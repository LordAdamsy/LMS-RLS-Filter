clc;clear;close all;

N = 500;
M = 11;
e_mean = zeros(1,N);
e_r_mean = zeros(1,N);
mu = 0.1; 

for i = 1:20
    x = 2*(randn(1, N) > 0) - 1; 
    d = [zeros(1, 6), x(1:N-6)];
    
    h = [0.3, 0.9, 0.3];
    x_c = filter(h, 1, x);
    
    x_n = x_c + 0.001 * randn(1, N);
    
    W = zeros(1,M); 
    W_r = zeros(1,M);
    x_i = zeros(1,M); 
    x_r = zeros(1,M); 
    y = zeros(1,N);
    y_r = zeros(1,N);
    e = zeros(1,N); 
    e_r = zeros(1,N); 
    delta = 0.1;
    P = 1/delta*eye(M);
    lambda = 0.95;
    
    for k = 1:N
        x_i = [x_n(k) x_i(1:M-1)]; 
        y(k) = W*x_i.';
        e(k) =d(k) - y(k);
        W = W + mu*e(k)*x_i; 
    end
    
    for k = 1:N
        x_r = [x_n(k), x_r(1:M-1)];
        Pk = P * x_r.';
        K = Pk/(lambda + x_r*Pk);
        s = d(k) - W_r*x_r.';
        W_r = W_r + K.'*conj(s);
        P =  (1/lambda)*P - (1/lambda)*K*x_r*P;
        e_r(k)=s;
        y_r(k) = W_r*x_r.';
    end

    e_mean = e_mean + e;
    e_r_mean = e_r_mean + e_r;
end

e_mean = e_mean/20;
e_r_mean = e_r_mean/20;

figure
subplot(411)
plot(d, 'LineWidth', 1.1);
hold on;
ylim([-2, 2]);
title('expected');
subplot(412)
plot(x_n, 'LineWidth', 1.1);
hold on;
ylim([-2, 2]);
title('noise signal');
subplot(413)
plot(y, 'LineWidth', 1.1);
hold on;
ylim([-2, 2]);
title('LMS');
grid on;
subplot(414)
plot(y_r, 'LineWidth', 1.1);
hold on;
ylim([-2, 2]);
title('RLS');
grid on;

figure
semilogy(e_mean.^2', 'LineWidth', 1.5);
hold on;
semilogy(e_r_mean.^2, 'LineWidth', 1.5);
hold on;
legend('LMS', 'RLS');
xlabel('iter');
ylabel('MSE');
grid on;