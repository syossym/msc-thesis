
K = -10:0.01:10;
r1 = 0.75;
F1 = 4 * r1^2 / (1 - r1^2)^2;
I1 = 1 ./ (1 + F1 * sin(K /2).^2);
figure(1)

plot(K,I1, 'r');

hold on
r2 = 0.85;
F2 = 4 * r2^2 / (1 - r2^2)^2;
I2 = 1 ./ (1 + F2 * sin(K /2).^2);

plot(K , I2 , 'g');
hold on
r3 = 0.99;
F3 = 4 * r3^2 / (1 - r3^2)^2;
I3 = 1 ./ (1 + F3 * sin(K /2).^2);

plot(K , I3);
xlabel('\delta');
ylabel('T'); 
%title('The Transmitted Intensity of a Fabry-Perot Interferometer');
legend(['R=' num2str(r1)], ['R=' num2str(r2)], ['R=' num2str(r3)]);

figure(2)
r11=.5:.005:.99;
FF = 3.14 * sqrt(r11) ./ (1 - r11);
plot(r11,FF)
xlabel('R');
ylabel('F');
%title('Reflectivity Finesse versus Mirorr Reflectivity R');
figure(3)
FFC = 4 * r11.^2 ./ (1 - r11.^2).^2;
plot(r11,FFC)
ylim([-20 800]) 
xlabel('Mirror Reflectivity R');
ylabel('Contrast Factor FFC');
title('Contast Factor versus Mirorr Reflectivity R');
grid on


