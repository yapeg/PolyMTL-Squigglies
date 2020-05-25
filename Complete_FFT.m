global Rc SR Uc Ay lambda
SR = 1.15; FR = 0.3; LR = 4;

fc = 1; Rc = 1; Uc = 2*pi*fc*Rc;
Vp = Uc/SR;

%% No Forcing Part (0)
tspan = 0:1e-3:1e3; Y_0 = [1 0 pi/2];
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5 1e-5 1e-5], 'refine', 5);
[t,Y0] = ode23(@noforcing, tspan, Y_0, options);

X10 = Y0(:,1).*cos(Y0(:,2)); X20 = X10 + (Uc/SR).*(t(end) - t);
Y10 = Y0(:,1).*sin(Y0(:,2)); Y20 = Y10;

n = round(length(Y0)*0.5):length(Y0); m = round(length(Y0)*0.9):length(Y0);
XF10 = X10(n); YF10 = Y10(n); figure(1); plot(XF10/Rc,YF10/Rc)
XF20 = X20(m); YF20 = Y20(m); figure(2); plot(XF20/Rc,YF20/Rc)
tF = t(n);

L = ceil(Uc*(tF(end) - tF(1))/Rc);
if mod(L,2) ~= 0
    L = L+1; end

% FFT on Y-coordinate
TFY0 = fft(YF10); P2Y0 = abs(TFY0);
P1Y0 = P2Y0(1:L/2+1); P1Y0(2:end-1) = 2*P1Y0(2:end-1); P1Y0 = P1Y0(2:end);
% FFT on X-coordinate
TFX0 = fft(XF10); P2X0 = abs(TFX0);
P1X0 = P2X0(1:L/2+1); P1X0(2:end-1) = 2*P1X0(2:end-1); P1X0 = P1X0(2:end);

f = 2*pi*(0:L/2)/L; f = f(2:end);
P1X0_n = P1X0/max([P1X0;P1Y0]);
P1Y0_n = P1Y0/max([P1X0;P1Y0]);

figure(3); plot(f,P1X0_n); hold on; plot(f,P1Y0_n); legend('wx','wy'); hold off

i = 1;
while P1Y0_n(i) < rms(P1Y0_n)
    i = i + 1; end
j = i;
while P1Y0_n(j) > rms(P1Y0_n)
    j = j + 1; end
[A,B] = max(P1Y0_n(i-1:j)); FR0 = f(i+B-2); fp = fc*FR0; 

%% Forcing part
fy = fp*FR; Ay = Rc*LR;
lambda = Uc/fy;

[t,Y] = ode23(@forcing, tspan, Y_0, options);
Yc = Ay*sin(2*pi*Uc*t/lambda);

X1 = Y(:,1).*cos(Y(:,2)); X2 = X1 + (Uc/SR).*(t(end) - t);
Y1 = Y(:,1).*sin(Y(:,2)) + Yc; Y2 = Y1;

XF1 = X1(n); YF1 = Y1(n); figure(4); plot(XF1/Rc,YF1/Rc)
XF2 = X2(m); YF2 = Y2(m); figure(5); plot(XF2/Rc,YF2/Rc)

% FFT on Y-coordinate
TFY = fft(YF1); P2Y = abs(TFY);
P1Y = P2Y(1:L/2+1); P1Y(2:end-1) = 2*P1Y(2:end-1); P1Y = P1Y(2:end);
% FFT on X-coordinate
TFX = fft(XF1); P2X = abs(TFX);
P1X = P2X(1:L/2+1); P1X(2:end-1) = 2*P1X(2:end-1); P1X = P1X(2:end);

P1X_n = P1X/max([P1X;P1Y]);
P1Y_n = P1Y/max([P1X;P1Y]);

figure(6); plot(f,P1X_n); hold on; plot(f,P1Y_n); legend('wx','wy'); hold off

% Extra
fun = @(t) Vp*sqrt(1+(2*pi*Ay*fy/Vp*cos(2*pi*fy*t)).^2);
lambda2 = integral(fun,0,1/fy);
if LR*FR*FR0 > pi/2
status = 'WARNING: Uc < Uw'; else       % Coiling speed lower than average wave speed
status = '-'; end
if lambda2 > lambda
status2 = 'WARNING: E < lambda'; else   % Extruded material within a wave lower than the own wave arc length
status2 = '-'; end

%%
function dY = noforcing(t,Y)
global Rc SR Uc
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = Uc*(cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)));
dY(2) = Uc*(sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2))) / Y(1);
dY(3) = Uc*kappa;
end
%%
function dY = forcing(t,Y)
global Rc SR Uc Ay lambda
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = Uc*(cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)) - 2*pi*Ay/lambda*cos(2*pi*Uc*t/lambda)*sin(Y(2)));
dY(2) = Uc*(sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2)) - 2*pi*Ay/lambda*cos(2*pi*Uc*t/lambda)*cos(Y(2))) / Y(1);
dY(3) = Uc*kappa;
end