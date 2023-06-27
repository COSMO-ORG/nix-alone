clear all;close all;clc;

load('./snow_data.txt')
load('./inp/WFJ_HS.txt')

ax_1 = subplot(8,1,1);
plot(snow_data(:,1),'r');hold on;plot(WFJ_HS,'g');
legend('Modeled snow height', 'measured snow height');  
hold off;

ax_2 = subplot(8,1,2);
plot(snow_data(:,2));
legend('Sensible heat flux');

ax_3 = subplot(8,1,3);
plot(snow_data(:,3));
legend('Latent heat flux');

ax_4 = subplot(8,1,4);
plot(snow_data(:,4));
legend('Albedo');

ax_5 = subplot(8,1,5);
plot(snow_data(:,5));
legend('air temperature');

ax_6 = subplot(8,1,6);
plot(snow_data(:,6));
legend('surface temperature');

ax_7 = subplot(8,1,7);
plot(snow_data(:,7));
legend('top layer thickness');

ax_8 = subplot(8,1,8);
plot(snow_data(:,8));
legend('water in the snow layer');

linkaxes([ax_1,ax_2,ax_3,ax_4,ax_5,ax_6,ax_7,ax_8],'x');
