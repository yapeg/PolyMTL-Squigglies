global Rc Ay SR LAMBDA

Tc = 1;           % steady coiling period                             (s)
Rc = 1;           % steady coiling radius                             (mm)
Uc = 2*pi/Tc*Rc;  % steady coiling speed                              (mm/s)

SR = 1.3;         % speed ratio  (Uc/Vp)
PR = 11;          % period ratio (Ty/Tc)
LR = 9;           % length ratio (Ay/Rc)

Vp = Uc/SR;       % platform speed                                    (mm/s)
Ty = Tc*PR;       % wave period                                       (s)
Ay = Rc*LR;       % wave amplitude                                    (mm)
                 
lambda = Uc*Ty;   % wave arc length                                   (mm)
fun = @(t) Vp*sqrt(1+(2*Ay*pi/Vp/Ty*cos(2*pi*t/Ty)).^2);
lambda2 = integral(fun,0,Ty);   % standard sine wave arc length       (mm)

S_int = 5*lambda; % total integration time                            (mm)

if lambda2 > lambda                   % the arc length of a perturbed flow
  output = 'Result might be invalid'  % is physically expected to be greater  
end                                   % than the one of a regular sine wave.

% ODE funtion resolution and ODE options
options = odeset('RelTol', 1e-5, 'AbsTol', [1e-5 1e-5 1e-5], 'refine', 5);
[S,Y] = ode23(@sine, [0 S_int], [1 0 pi/2], options);

% Sinusoidal forcing
Yc = Ay*sin(2*pi*S/LAMBDA);

% Contact point trace
X1 = Y(:,1).*cos(Y(:,2));
Y1 = Y(:,1).*sin(Y(:,2)) + Yc;

% Deposited trace
X2 = X1 + Vp/Uc.*(S(length(S)) - S);
Y2 = Y1;


%% Plotting (optional)
% Text box
txt = cell(3,1); txt{1} = ['SR = ', num2str(SR)]; txt{2} = ['PR = ', num2str(PR)];
                 txt{3} = ['LR = ', num2str(LR)];

% Creation of a folder containing this and other plots done, named after its dynamic characteristics
fname = ['SR', num2str(SR), 'PR', num2str(PR), 'LR', num2str(LR)];
mkdir('Plots',fname)
fpath = fullfile(pwd,'Plots',fname);

% Plot: contact point trace
figure(1)
plot(X1,Y1)
title('Sinusoidal movement - Plot of contact point trace')
xlabel('X(s) (mm)', 'FontWeight', 'Bold')
ylabel('Y(s) (mm)', 'FontWeight', 'Bold')
axis([min(X1) max(X1), min(Y1) max(Y1)])
text(min(X1)+0.04*(max(X1)-min(X1)), max(Y1)-0.10*(max(Y1)-min(Y1)), txt, 'FontWeight', 'Bold', 'FontSize', 9)
grid
% Saves figure in the mentioned folder
saveas(figure(1),fullfile(fpath,'plot1.fig'))

% Plot: deposited trace
figure(2)
plot(X2,Y2)
title('Sinusoidal movement - Plot of deposited trace')
xlabel('X(s) (mm)', 'FontWeight', 'Bold')
ylabel('Y(s) (mm)', 'FontWeight', 'Bold')
axis([min(X2) max(X2), min(Y2) max(Y2)])
text(min(X2)+0.04*(max(X2)-min(X2)), max(Y2)-0.10*(max(Y2)-min(Y2)), txt, 'FontWeight', 'Bold', 'FontSize', 9)
grid
% Saves figure in the mentioned folder
saveas(figure(2),fullfile(fpath,'plot2.fig'))


%% Function
% Let Y(1) = r, Y(2) = psi, Y(3) = theta, and dY its derivatives with respect to S (arc length along the deposited trace)
function dY = sine(S,Y)
global Ay Rc SR LAMBDA
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)) - 2*pi*Ay/LAMBDA*cos(2*pi*S/LAMBDA)*sin(Y(2));
dY(2) = (sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2)) - 2*pi*Ay/LAMBDA*cos(2*pi*S/LAMBDA)*cos(Y(2))) / Y(1);
dY(3) = kappa;
end
