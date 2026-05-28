function [exportData] = Wave_Forces(depth,height,period,diameter,z_c)

%% Inputs
%
%   depth - water depth in the wave region (m)
%   
%   height - wave height (m)
%
%   period - wave period (s)
%
%   diameter - cylinder diameter (m)
%
%   z_c - height on the cylinder to analyse forces at (m)

%% Initialize know variables
g = 9.807;
rho = 1025;
%z_c = 5;
%z_c = -depth;


%% Validate input data

    % Check if inputs are positive values
    if any([depth, height, period, diameter] <= 0)
        error('All inputs must be positive values.');
    end
    
    % Check if inputs are numerical
    if (isa([depth, height, period, diameter],"double") == 0)
        error('Inputs must be numerical values.');
    end

% Calculate angular frequency (omega)
omega = 2 .* pi ./ period;


% Solve for wave number (k) using angular frequency and make positive
% Make a guess for value of k using wave_length function
waveLength = Wave_Length(period,depth);
k = zeros(1,length(period));
for i = 1:length(period)
    k(i) = 2.*pi./waveLength(i); %kguess
    
   % dipf=@(k) omega(i).^2 - g.*k.*tanh(k.*depth);
    
   % k(i) = fsolve(dipf,kguess);
end
%k = abs(k);

% Calculate max velocity (vm)
vm = ((g.*height)./(2.*omega)).*k;

% Calculate Keulegan-Carpenter (kc) number
kc = (vm.*period)./diameter;

% Calculate reynolds number (re)
v = 1.04.*(10.^(-6));

re = (vm.*diameter)./v;

% Parse wave data to Wave_Theory function to get correct theory
%[theory] = Wave_Theory(height,depth,period);


%% Load inertia coefficients graph and set scale then plot the kc number as a
% vertical line
figure
hold on
inertia = imread("Inertia_CoEf.png");
imagesc([0 50],[1.4 2.1],inertia,AlphaData=0.9);
title('Inertia Coefficient Vs Keulgan-Carpenter Number');
xlabel('KC');
ylabel('Cm');

xlim([0 50])
ylim([1.4 2.1])
% request user input to select correct coefficient.
C_m = zeros(1,length(kc));
for i = 1:length(kc)
    xline(kc(i),'r','LineWidth',1.5)
    [~,C_m(i)] = ginput(1);
end
scatter(kc,C_m,'diamond','filled','b')
hold off

% Repeat for drag coefficients
figure
hold on
drag = imread("Drag_CoEf.png");
imagesc([0 50],[0.5 1.6],drag,AlphaData=0.9);
title('Drag Coefficient Vs Keulgan-Carpenter Number');
xlabel('KC');
ylabel('Cd');

xlim([0 50])
ylim([0.5 1.6])
% Request user input for correct coeff
C_d = zeros(1,length(kc));
for i = 1:length(kc)
    xline(kc(i),'r','linewidth',1.5)
    [~,C_d(i)] = ginput(1);
end
scatter(kc,C_d,'diamond','filled','b','linewidth',2)
hold off

%% Calculate the required values to plot wave on the different forces
% regimes graph to determine inertia or drag dominance
HD = height./diameter;
DF_XAxis = (pi.*diameter)./waveLength;

figure
hold on
DiffForces = imread("Different _Forces.png");
title('Different Forces Regimes');
xlabel('\piD/\lambda');
ylabel('H/D');

xscale("log");
yscale("log");
imagesc([0.01 10],[0.01 100],DiffForces,AlphaData=0.9);

scatter(DF_XAxis,HD,'filled','r')
% x axis 0.48 for difraction/inertia dominance
ylim([0.01 100])
xlim([0.01 10])
axis('square')

hold off

%% Determin if wave has inertia or drag dominance and calculate force

for i = 1:length(DF_XAxis)
    if DF_XAxis(i) > 0.48
    
        % Defraction dominance, use MacCamy-Fuchs equation
    
        x(i) = k(i).*(diameter/2);
    
        for h = 0:4
            
            J(h+1) = besselj(h,x(i));
            Y(h+1) = bessely(h,x(i));
        
        end
    
        Jd_1 = (1/x(i)) .* J(2)-J(3);
        Yd_1 = (1/x(i)) .* Y(2) - Y(3);
    
        alpha(i) = atan(Jd_1/Yd_1);
    
        A(i) = (J(2).^2 + Y(2).^2).^(-1/2);
    
        max_inertia_F(i) = -2 .* ((rho.*g.*height(i))./k(i).^2) .* A(i) .* tanh(k(i).*depth);
        max_Drag_F(i) = 0;
    
        max_inertia_M(i) = 2 .* ((rho.*g.*height(i))./k(i).^3) .* A(i) .* ((k(i).*z_c.*sinh(k(i).*depth)+cosh(k(i).*depth)-1)/cosh(k(i).*depth));
        max_Drag_M(i) = 0;

        x1 = 0;
        x2 = period(i);

        F = @(t) -2 .* ((rho.*g.*height(i))./k(i).^2) .* A(i) .* tanh(k(i).*depth) .* tanh(omega.*t - alpha);
        M = @(t) 2 .* ((rho.*g.*height(i))./k(i).^3) .* A(i) .* ((k(i).*z_c.*sinh(k(i).*depth)+cosh(k(i).*depth)-1)/cosh(k(i).*depth)) .* sin(omega.*t - alpha);

        [~,Fm] = fminbnd(F,x1,x2);
        [~,Mm] = fminbnd(M,x1,x2);

        FMax(i) = abs(Fm);
        MMax(i) = abs(Mm);
    
    elseif DF_XAxis(i) < 0.48
    
        % Inertia Dominance, uses morrison forces
        % Calculate maxima
   
        max_inertia_F(i) = -C_m(i) .* ((pi.*rho.*(diameter.^2).*(omega(i).^2).*height(i))./(8.*k(i)));
        
        max_Drag_F(i) = C_d(i) .* ((rho.*diameter.*omega(i).^2.*height(i).^2)./(16.*(sinh(k(i).*depth)).^2)) .* ((sinh(2.*k(i).*depth)./(2.*k(i)))+depth);
        
        F = @(t) (-C_m(i) .* ((pi.*rho.*(diameter.^2).*(omega(i).^2).*height(i))./(8.*k(i))) .* sin(omega(i).*t)) + (C_d(i) .* ((rho.*diameter.*(omega(i).^2).*(height(i).^2))./(16.*((sinh(k(i).*depth)).^2))) .* ((sinh(2.*k(i).*depth)./(2.*k(i)))+depth) .* cos(depth.*t) .* abs(cos(depth.*t)));
        
        max_inertia_M(i) = C_m(i) .* ((pi.*rho.*(diameter.^2).*(omega(i).^2).*height(i))./8) .* (((2.*(sinh((k(i).*depth)./2).^2))./((k(i).^2).*sinh(k(i).*depth)))+(z_c./k(i)));
        
        max_Drag_M(i) = C_d(i) .* ((rho.*diameter.*(omega(i).^2).*(height(i).^2))./(16.*(sinh(k(i).*depth).^2))) .* (((1-cosh(2.*k(i).*depth))./(4.*(k(i).^2))) - (z_c.*((sinh(2.*k(i).*depth))./(2.*k(i)))) - (depth.*z_c) - (depth.^2)./2);
        
        M = @(t) (C_m(i) .* ((pi.*rho.*(diameter.^2).*(omega(i).^2).*height(i))./8) .* (((2.*(sinh((k(i).*depth)./2).^2))./((k(i).^2).*sinh(k(i).*depth)))+(z_c./k(i))) .* sin(omega(i).*t)) + (C_d(i) .* ((rho.*diameter.*(omega(i).^2).*(height(i).^2))./(16.*(sinh(k(i).*depth).^2))) .* ((1-cosh(2.*k(i).*depth)./(4.*(k(i).^2))) - (z_c.*((sinh(2.*k(i).*depth))./(2.*k(i)) - (depth.*z_c) - (depth.^2)./2))) .* cos(depth.*t) .* abs(cos(depth.*t)));
        
        x1 = 0;
        x2 = period(i);
        
        [~,Fm] = fminbnd(F,x1,x2);
        [~,Mm] = fminbnd(M,x1,x2);
        
        FMax(i) = abs(Fm);
        MMax(i) = abs(Mm);
    
    end
end


%% Export and display data 

% Export key values in tabular format
exportData = table(height',k',vm',kc',re',C_m',C_d',max_Drag_F',max_inertia_F',FMax',max_Drag_M',max_inertia_M',MMax','VariableNames',["Wave Height","Wave number","Orbital Velocity","Keulegan-Carpenter Number","Reynolds Number","Inerta Coefficient","Drag Coefficient","Maximum Drag Force (N)","Maximum Inertia Force (N)","Maximum Force (N)","Maximum Drag Moment (Nm)","Maximum Inertia Moment (Nm)","Maximum Over-Turning Moment (Nm)"]);

% Display as a figure
figure
% Get the table in string form.
wheString = evalc('disp(exportData)');
% Use TeX Markup for bold formatting and underscores.
wheString = strrep(wheString,'<strong>','\bf');
wheString = strrep(wheString,'</strong>','\rm');
wheString = strrep(wheString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',wheString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Write export data to file

name = sprintf("WaveForces_%s.dat",datetime('now'));
writetable(exportData,name);

end