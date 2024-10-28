function [f_l, f_h] = PlotWavefunction(E,k_t,p,c);

% Global variables
global gamma1  gamma2  gamma3  V  N  dz z;

tf_total = eye(4);
tf       = zeros(4,4);
tf(1,3)  = 1;
tf(2,4)  = 1;

f_h(1)   = 0;
f_l(1)   = 0;
f_h(2)   = p;
f_l(2)   = c;

f = [0, 0, p, c]';

% figure(3);
% plot(f_final(3,:));
% hold on;
% plot(f_final(4,:));

for (n=2:N-1)

    if (gamma1(n)~=gamma1(n-1) | n==2)

        C1       = gamma1(n)^2 - 4*gamma2(n)^2;
        Cf       = 1 + 3*(k_t^2)*(dz^2)*(gamma3(n)^2)/C1;

        tf(3,1)  = (-1 + 3*(k_t^2)*(dz^2)*(gamma3(n)^2)/C1) / Cf;
        tf(3,2)  = (2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)-2*gamma2(n))) / Cf;

        tf_33    = ( 2 + (k_t^2)*(dz^2)*(gamma1(n)+gamma2(n))/(gamma1(n)-2*gamma2(n)) - 3*(k_t^3)*(dz^3)*gamma2(n)*gamma3(n)/C1 ) / Cf;
        tf_34    = ( sqrt(3)*(k_t^2)*(dz^2)*gamma2(n)/(gamma1(n)-2*gamma2(n)) - 2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)-2*gamma2(n)) - sqrt(3)*(k_t^3)*(dz^3)*gamma3(n)*(gamma1(n)-gamma2(n))/C1 ) / Cf;

        tf(4,1)  = (-2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)+2*gamma2(n))) / Cf;
        tf(4,2)  = (-1 + 3*(k_t^2)*(dz^2)*(gamma3(n)^2)/C1) / Cf;

        tf_43    = ( sqrt(3)*(k_t^2)*(dz^2)*gamma2(n)/(gamma1(n)+2*gamma2(n)) + 2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)+2*gamma2(n)) + sqrt(3)*(k_t^3)*(dz^3)*gamma3(n)*(gamma1(n)+gamma2(n))/C1) / Cf;
        tf_44    = ( 2 + (k_t^2)*(dz^2)*(gamma1(n)-gamma2(n))/(gamma1(n)+2*gamma2(n)) + 3*(k_t^3)*(dz^3)*gamma2(n)*gamma3(n)/C1) / Cf;

        gammap   = gamma1(n) + 2*gamma2(n);
        gammam   = gamma1(n) - 2*gamma2(n);
        const34  = sqrt(3)*k_t*(dz^3)*gamma3(n)/C1;

    end

    dE       = V(n) - E;

    tf(3,3)  = tf_33 + ( (dz^2)*dE/gammam ) / Cf;
    tf(3,4)  = tf_34 + ( -const34*dE ) / Cf;
    tf(4,3)  = tf_43 + ( const34*dE ) / Cf;
    tf(4,4)  = tf_44 + ( (dz^2)*dE/gammap ) / Cf;

    tf;
    tf_total = tf * tf_total;
    f        = tf_total*[0,0,p,c]';
    %f = tf * f;
    
    f_h(n)   = f(1);
    f_l(n)   = f(2);
    f_h(n+1) = f(3);
    f_l(n+1) = f(4);

%     figure(3);
%     drawnow;
%     plot(f_h);
%     hold on;
%     plot(f_l, 'r');
%     title(['E=' num2str(E)]);

    if (n==N-10)
        tf_p = tf_total;
    end
end

figure(3);
%plot(z, f_l, 'b', z, f_h, 'r');
%hold on;
 
norm    = sqrt(trapz(z, f_l.^2+f_h.^2));

f_l     = f_l./norm;
f_h     = f_h./norm;
plot(z, f_l, 'g', z, f_h, 'y');
hold off;


