function [xx] = resampling(type,N,waga,x_p)
switch type
    case 'systematic'
        j=1; sumQ=waga(j);
            u=rand()/N;
            for i=1:N
                while sumQ < u
                    j=j+1;
                    sumQ=sumQ+waga(j);
                end
                xx(:,i)=x_p(:,j);
                u=u+1/N;
            end
    case 'stratified'
    j=1; sumQ=waga(j);
    for i=1:N
        u=(rand()+i-1)/N;
        while sumQ<u
            j=j+1;
            sumQ=sumQ+waga(j);
        end
        xx(:,i)=x_p(:,j);
    end
    case 'residual'
    Nr=N; jj=0;
        for i=1:N
            j=floor(waga(i)*N);
            for ii=1:j
                jj=jj+1;
                xx(:,jj)=x_p(:,i);
            end
            waga(i)=waga(i)-j/N;
            Nr=Nr-j;
        end
        if Nr>0
            for i=1:N
                waga(i)=waga(i)*N/Nr;
            end
            j=1; sumQ=waga(1);
            for i=1:Nr
                u=(rand()+i-1)/Nr;
                while sumQ<u
                    j=j+1;
                    sumQ=sumQ+waga(j);
                end
                jj=jj+1;
                xx(:,jj)=x_p(:,j);
            end
        end
    otherwise
        error('unexpected resmpling type, aborting')
    end
end