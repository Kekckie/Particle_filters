function [aRMSE, Jx, Jy, eps] = run_particle_filter(M, Vd, Vn, N, f, h, x_pocz, resampling_type)
% run filter
% Returns values of indicators given object and its parameters, number of
% iterations and number of particles used for prediciton

% step_1: Inicialization
    x = x_pocz;
    x_hat = x_pocz; %initial values of estimated state variables
    for t=2:M %state variable init
        x(:,t) = f(x(:,t-1))' + Vd .* randn(size(x,1),1);
    end
    for t=1:M %output values init
        y_true(:,t) = h(x(:,t))';
        y_pomiar(:,t) = h(x(:,t))' + Vn .* randn(size(y_true,1),1);
    end
    for p=1:size(N,2)
        for i=1:N(p) % particle init.
            x_p(:,i) = [x_pocz];
            waga(i) = 1/N(p);
        end
        for t=2:M % step_2: prediction
            for i=1:N(p)
                x_p_upd(:,i) = f(x_p(:,i))' + Vd .* randn(size(x,1),1);
            end
            x_p = x_p_upd;

            for i=1:N(p) % step_3: actualization
                y_obl = h(x_p(:,i))';
                rozn = y_pomiar(:,t) - y_obl;
                waga(i) = (1/sqrt(2*pi*Vn(1)^2)) * exp(-0.5*(rozn(1)^2/(Vn(1)^2)));
                for j=2:size(y_true,1)
                    waga(i) = waga(i) * (1/sqrt(2*pi*Vn(j)^2)) * exp(-0.5*(rozn(j)^2/(Vn(j)^2)));
                end     
            end

           %draw_particles() %<------------------drawin particles and their weights

            suma_QQ = sum(waga); % step_4: normalization
            for i=1:N(p)
                waga(i) = waga(i)/suma_QQ;
            end
            % step_5 resampling method 12 systematic
            
            x_p=resampling(resampling_type,N(p),waga,x_p);
            for i=1:N(p)
                waga(i) = 1/N(p);
            end

            % step_6: estimation of state variables
            for i=1:size(x_p,1)
                x_hat(i,t) = mean(x_p(i,:));
            end
       

        %draw_particles() %<------------------drawin particles and their weights

        end
        for i=1:M %estimation of outputs based on estimation of state variables
                y_hat(:,i) = h(x_hat(:,i))';
        end
        for i=1:size(x,1)
            aRMSEj(i) = sqrt((1/M)*(sum((x_hat(i,:)-x(i,:)).^2)));
        end
        aRMSE(p) = (1/size(x,1))*(sum(aRMSEj));
        % Relative Mean Square Error Jx
        Jx(p) = 0;
        for i=1:size(x,1)
            Jxj(i) = (1/(M*(Vd(i)^2))*(sum((x_hat(i,:)-x(i,:)).^2)));
        end
        Jx(p) = (1/size(x,1))*(sum(Jxj));
        % Relative Mean Square Error of measurments Jy
        Jy(p) = 0;
        for i=1:size(y_true,1)
            Jyj(i) = (1/(M*(Vn(i)^2)))*(sum((y_hat(i,:)-y_true(i,:)).^2));
        end
        Jy(p) = (1/size(y_true,1))*sum(Jyj);
        % Relation of estimation error to the measure error epsilon
        eps(p) = 0;
        top = 0;
        bot = 0;
        for i=1:size(y_true,1)
            for j=1:M
                top = sum(abs((y_hat(i,:)-y_true(i,:))));
                bot = sum(abs((y_pomiar(i,:)-y_true(i,:))));
            end
        end
        eps(p) = (1/size(y_true,1))*(top/bot);
    end
    % plot(1:M,y_hat(1,:))
    % hold on
    % plot(1:M,y_true(1,:))
end

