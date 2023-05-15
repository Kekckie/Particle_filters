clear, clc, close all
M = 100;
Vd = 0.5;
Vn = 1;
N = 100;
f=@(x)[x(1),x(2)];
h=@(x)[x(1)+x(2),x(1)-x(2),x(1)*x(2)];
x = [0.1;0.1];
x_hat = x;
%inicjalizacja_cząsteczek==================================================
for i=1:N
    x_p(:,i) = [x];
    waga(i) = 1/N;
    waga_prev = waga;
end
for t=2:M
    x(:,t) = f(x(:,t-1))' + Vd * randn(2,1);
end
for t=1:M
    y_true(:,t) = h(x(:,t))';
    y_pomiar(:,t) = h(x(:,t))' + Vn * randn(3,1);
end
for t=2:M
    %krok_2_predykcja===========================
    for i=1:N
        x_p_upd(:,i) = f(x_hat(:,t-1))' + Vd * randn(2,1);
    end
    x_p = x_p_upd;
    
    %krok_3_aktualizacja========================
    for i=1:N
        y_obl = h(x_p(:,i))';
        rozn = y_pomiar(:,t) - y_obl;
        wagai(i) = ((1/sqrt(2*pi*Vn*Vn)) * exp(-0.5*(rozn(1)^2/(Vn*Vn))) ...
            * (1/sqrt(2*pi*Vn*Vn)) * exp(-0.5*(rozn(2)^2/(Vn*Vn))) ...
            * (1/sqrt(2*pi*Vn*Vn)) * exp(-0.5*(rozn(3)^2/(Vn*Vn))))*waga_prev(i);      
    end
    %krok_4_normalizacja
    suma_QQ = sum(wagai);
    for i=1:N
        wagai(i) = wagai(i)/suma_QQ;
    end
%==========RESAMPLING_PROBA_2=======================================
%     index = randi(N-1);
%     beta = 0;
%     mw = max(x_p_upd(3,:));
%     
%     for i=1:N
%        beta = beta + 2*mw*rand();
%        while beta > x_p_upd(3,i)
%           beta = beta - x_p_upd(3,i);
%           index = rem((index + 1),N);
%           if index == 0
%               index = 1;
%           end
%        end
%        new_p(:,i) = x_p_upd(:,index);
%     end
%     x_p = new_p;
%     x_p(3,:) = 1/N;
%==========================RESAMPLING_102_alg2==================
%     c(1) = 0;
%     for i=2:N
%         c(i) = c(i-1) + waga(i);
%     end
%     i = 1;
%     c(100) = 1;
%     u(1) = rand()/N;
%     for j=1:N
%        u(j) = u(1) + (1/N)*(j-1);
%        while u(j) > c(i)
%            i = i + 1;
%        end
%        x_p(:,j) = x_p_upd(:,i);
%        waga(i) = 1/N;
%     end

    % resampling-12===========================
    j=1; sumQ=wagai(j);
    u=rand()/N;
    for i=1:N
        while sumQ < u
            j=j+1;
            sumQ=sumQ+wagai(j);
        end
        xx(:,i)=x_p(:,j);
        u=u+1/N;
    end
    x_p=xx;
    for i=1:N
        wagai(i) = 1/N;
    end
    for j=1:N
        % dobranie xj i rzeczy do obliczenia wagi xj, tak jak w SIR
        y_obl = h(x_p(:,i))';
        rozn = y_pomiar(:,t) - y_obl;
        xj(:,j) = f(x_p(:,j))' + Vd * randn(2,1);
        wagaj(j) = ((1/sqrt(2*pi*Vn*Vn)) * exp(-0.5*(rozn(1)^2/(Vn*Vn))) ...
            * (1/sqrt(2*pi*Vn*Vn)) * exp(-0.5*(rozn(2)^2/(Vn*Vn))) ...
            * (1/sqrt(2*pi*Vn*Vn)) * exp(-0.5*(rozn(3)^2/(Vn*Vn))));
        % normalizacja wagi xj
        suma_QQ = sum(wagai);
        for i=1:N
            wagai(i) = wagai(i)/suma_QQ;
        end
        % nowa waga liczona jako stosunek między wagą tą wyżej a niżej
            waga(j) = wagai(j)/wagaj(j);
    end
    waga = waga_prev;
    %krok_6_estymacja========================
    % draw_particles();
    x_hat(1,t) = mean(x_p(1,:));
    x_hat(2,t) = mean(x_p(2,:)); 
    

end
% plot(1:N,x_hat(1,:))
% hold on
% plot(1:N,x_hat(2,:))

plot(1:M,x_hat(2,:))
hold on
plot(1:M,x(2,:))
legend('estm','true')
% plot(1:N,y(1,:))



