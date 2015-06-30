% Saccade force-velocity estimates. This is a rough planar estimate assuming 
% two linear actuators at distance r from centre of rotation; centre of 
% rotation at centre of mass. The rods are assumed to be long, so their
% angle doesn't change with the camera angle. 

dt = .001;
time = 0:dt:.2;

% example saccade using form from Opstal & Gisbergen (1987) and reasonable
% large parameters ... 
beta = .015;
gamma = 5;
omega = 130 * (time/beta).^(gamma-1) .* exp(-time/beta) * pi/180; % angular velocity
alpha = [0 diff(omega)/dt]; % angular acceleration
relTheta = cumsum(dt*omega); %angle relative to starting angle

% guesses for camera-lens assembly ... 
mass = .083 + .078 + 2*.035;
I = mass * .05^2;

figure, set(gcf, 'Position', [312   354   865   420])

subplot(1,2,1), set(gca, 'FontSize', 18)
plot(time, omega*180/pi)
xlabel('Time (s)', 'FontSize', 18)
ylabel('Angular Velocity (deg/s)', 'FontSize', 18)

subplot(1,2,2), set(gca, 'FontSize', 18)
hold on
c = {'r', 'b', 'g', 'c'};

plot(0:15, .25-.25/15*(0:15), 'k') % see http://www.pi-usa.us/products/PDF_Data/U264_High_Speed_Piezo_Motor_Actuator.pdf
xlabel('Force (N)', 'FontSize', 18)
ylabel('Velocity (m/s)', 'FontSize', 18)

plot(0,0,c{1}), plot(0,0,c{2}), plot(0,0,c{3}), plot(0,0,c{4}); % dummy points for legend 
legend('Max', '.005m', '.01m', '.015m', '.02m')

r = [.005 .01 .015 .02]; %radii (note: we also get better angular resolution with longer r)
theta0 = [-pi/6 0 pi/6]; %camera angle at start of saccade (0 is straight forward)
for i = 1:length(r)
    for j = 1:length(theta0)
        theta = theta0(j) + relTheta;
        fp = I*alpha/r(i); %perpendicular force on camera needed to produce acceleration alpha
        f = fp./cos(theta); %corresponding force applied by actuator
        v = r(i).*cos(theta).*omega; %linear velocity of actuator 
        plot(f/2, v, c{i}) %divide by two to plot force for each actuator
    end
end
