function  [] = MaxwellBoltz(vel, tset)

%This function can be used to verify if a proper Maxwell-Boltzmann
%distribution is achieved

kb = 1.3806485279e-23; %J/K
amass = 26.9815385*1.66053904020e-27; %kg

vs = [0:0.1:1000];

velxys = [-1000:0.1:1000];

%3D!
%max_boltz = 4*pi*vs.^2*(amass/(2*pi*kb*tset))^(3/2).*exp(-amass*vs.^2/(2*kb*tset));

%2D
max_boltz =amass*vs/(kb*tset).*exp(-amass*vs.^2/(2*kb*tset));

sigma = sqrt(kb*tset/amass);
mu=0;
normal_distribution_directions = 1/(sqrt(2*pi*sigma^2))*exp(-(velxys-mu).^2/(2*sigma^2));

vels = sqrt(vel(:,1).^2 + vel(:,2).^2);

figure(20)
subplot(2,2,1)
histogram(vels,50,'Normalization','pdf')
hold on
plot(vs, max_boltz)
xlim([vs(1), vs(end)])
ylim([min(max_boltz), 1.25*max(max_boltz)])
ylabel('Probability')
xlabel('Speed')
hold off

subplot(2,2,3)
histogram(vel(:,1),50,'Normalization','pdf')
hold on
plot(velxys, normal_distribution_directions)
xlim([velxys(1), velxys(end)])
ylim([min(normal_distribution_directions), 1.25*max(normal_distribution_directions)])
hold off
xlabel('Velocity in x-direction')
ylabel('Probability')

subplot(2,2,4)
histogram(vel(:,2),50,'Normalization','pdf')
hold on
plot(velxys, normal_distribution_directions)
xlim([velxys(1), velxys(end)])
ylim([min(normal_distribution_directions), 1.25*max(normal_distribution_directions)])
hold off
xlabel('Velocity in y-direction')
ylabel('Probability')

drawnow
end

