clc; clear; close all;
%% User guide
%titles of the excel files for each burn will have to be
%"propellant_test_data_1","propellant_test_data_2" "propellant_test_data_3" 
%or the name can be changed in Load Data.
%until update, length of P and t data of each file will need to be the same
%script is currently set up for both ends uninhibited, tubular geometry.
%this can be changed by modifying Ab_0 and the %moooove section of 
%findRegression
%for to get a and n into steady region, select steady state pressure, P_ss,
%from the graphs.


%% constants
%Motor Geometry, currently test data
Do_0 = 43.1;                %mm      initial outer diameter
Di_0 = 13.88;               %mm      initial inner diameter 
L0   = 65;                  %mm      initial length
Ro_0 = Do_0/2;              %mm      initial outer radius
Ri_0 = Di_0/2;
% Ro_0 = ;                  %If grain is measured in radius
% Ri_0 = ;

numGrains = 2; %Number of grains per firing, so long as the total length 
%of the propellant is correct the formulas should work. 
numFires = 3;

%Test Nozzle Geometry in Imperial
D_tImp = [5/16,21/64,11/32,23/64];
D_tMet = D_tImp*25.4; %immigrate to Mother Russia's glorious communist units, %25.4mm to 1 in
A_t = (pi/4)*((D_tMet(:)).^2);
A_t(1) = 81.07; %for test data for first data set to ensure working on % other computer.
rho  = 1.912;               %g/cm^3 density of propellant
m_p  = 0.324;               %kg      mass of propellant


Ab_0 = numGrains*(2*pi*Ri_0*L0 + 2*pi*(Ro_0^2 -Ri_0^2)); %inititialize burning area for both ends uninhibited
% Ab_0 = numGrains*(2*pi*Ri_0*L0 + 1*pi*(Ro_0^2 -Ri_0^2)); %inititialize burning area for one ends uninhibited
% Ab_0 = numGrains*(2*pi*Ri_0*L0); inhibited ends
% Kn0(1,numFires)= Ab_0/A_t; 

%% Load Data (accepts .xlsx or .csv)
%loops through each file, pulling the data and filling them into pressure
%and time matrix
for k = 1:numFires
    filename = sprintf("propellant_test_data_%d",k);
    filenameExcel = append(filename,".xlsx");
    filenameCSV = append(filename,".csv");

    if exist(filenameExcel, 'file')
        data = readtable(filenameExcel,'PreserveVariableNames',true);
    elseif exist(filenameCSV,'file')
        data = readtable(filenameCSV,'PreserveVariableNames',true);
    else
        error('No data file found. Put propellant_test_data_%d.xlsx or .csv in the current folder.',k);
    end
    % Time column (expects "t (s)" or "time")
    if any(strcmp(data.Properties.VariableNames,'t_(s)'))
        t_k = data.('t_(s)');
    elseif any(strcmp(data.Properties.VariableNames,'time'))
        t_k = data.time;
    else
        error('Time column not found. Expected header "t (s)" or "time".');
    end
    k;
    t(:,k) = t_k;


    % Pressure column (prefer SI)
    if any(strcmp(data.Properties.VariableNames,'Po_SI_(N/m^2)'))
        P_k = data.('Po_SI_(N/m^2)');
    elseif any(strcmp(data.Properties.VariableNames,'Po_(psi)'))
        P_k = data.('Po_(psi)') * 6894.75729; % psi -> Pa
    else
        error('Pressure column not found. Expected "Po_SI (N/m^2)" or "Po (psi)".');
    end
    k;
    P(:,k) = P_k(:);
end
N = length(t(:,numFires)); %finds no. of datapoints

%% run functions over the four sets of data process on P and t data, 
% returns r[mm/s] and other values for checking regression
  c = zeros(N,numFires);
  c = getCharacteristicVelocity(P,t,A_t,m_p,numFires);
  [Ri,Ro,L,Ab,ds,r,Kn] = findRegression(P,t,rho,A_t,Ro_0,Ri_0,L0,Ab_0,numGrains,c,numFires,N);


%% Powerfit
% requires curve fitting toolbox
% select pressure a for steady state applies power fit to r vs P to get a % and n 
%Graph for selecting
    figure('Name','Pressure Vs time'); hold on;
    for i=1:numFires
    plot(t, P);
    xlabel('Time (s)');
    ylabel('Pressure (Pa)');
    title('Pressure vs Time for all fires, find steady state pressures');
    legendTitles{i} = sprintf("Static Fire %d",i);
    end
    legend(legendTitles);
 

% P_ss = [3.1,5.8,7.6]*10^6
P_ss = [3.1,3.1,3.1,3.1]*10^6; %steady state region defined by P_ss
% Adjust the steady state pressures for each firing from pressure graph
%for the individual fires
  sprintf("Now printing data for each static fire")
for l =1:numFires
    l;
    r_l = r(:,l); %r has to be a column vector for fit function
    %takes indices of first and last pressures above a minimum to get
    %steady state operation period
    idxLow = find(P(:,l) >= P_ss(l),1,'first'); 
    idxHigh = find(P(:,l) >= P_ss(l),1,'last');
    %fits the data to a r = aP^n with powerfit.a as a and powerfit.b = n
   
        sprintf("Power fit for Static Fire %d",l)
    [powerfit] = fit(P(idxLow:idxHigh,l)*10^-6, r_l(idxLow:idxHigh),'power1')
    P_linear = P(idxLow:idxHigh,l); %the pressure region which is above steady state Pressure, and is in linear region which Saint Roberts Law
    r_linear = r_l(idxLow:idxHigh);
 

    burn_coefficient(l) = powerfit.a;
    pressure_exponent(l) = powerfit.b;
    confidenceA = confint(powerfit,0.95);
    r_StRobLaw = burn_coefficient(l)*((P_linear).^(pressure_exponent(l))) %takes r = aP^n fit to predict regression from pressure
    r_StRobLawMPa = burn_coefficient(l)*((P_linear*(10^-6)).^(pressure_exponent(l)))
    figure
    periwinkle = '#CCCCFF'; %The lion cares meticulously of the color of the line
    loglog(P_linear,r_linear,'bo'); hold on; %plots points
    loglog(P_linear, r_linear, 'r-', 'LineWidth', 1, 'Color',periwinkle); %plots line between points
    % loglog(xAxis, burn_coefficient(l)*(xAxis).^(pressure_exponent(l)), 'bo');
    loglog(P_linear, r_StRobLawMPa, 'r-', 'LineWidth', 1, 'Color','Magenta');
    title(sprintf(['St. Robert''s Law of static fire %d Fit(log–log): n = %.3f,' ...
        ' a = %.3g [mm/s]'], l, pressure_exponent(l), burn_coefficient(l)));
    subtitle(sprintf(['95 Percent CI, a = (%.3f,%.3f) [mm/s], n = (%.3f,%.3f)'],confidenceA(1),confidenceA(2),confidenceA(3),confidenceA(4)))
    hold off
    
end


%% functions
function [c] = getCharacteristicVelocity(P,t,A_t,m_p,numFires) %c is velocity at end of motor
   for k = 1:numFires
    %equation: c = (A_throat/mass_propellant)*ΣPressure*Δt
    P_sum = sum(P(:,k));  %pa
    dt = mean(diff(t(:,1))); %s
    A_t(k);
    c(k) = ((A_t(k)*(10^-6))/m_p)*P_sum*dt; %The test data was modified σ_σ
    end
end

%P,t,A_t are going to be passed into 

function [Ri,Ro,L,Ab,ds,r,Kn] = findRegression(P,t,rho,A_t,Ro_0,Ri_0,L0,Ab_0,numGrains,c,numFires,N) %inputs, (p,t), constants 

%Iterative solving loop because Ab changes with time as well as Δs 
%fzero to find Δs - P(t)/Kn*Rho*c_star = 0, value to find is Δs then r
%initialize Ro, Ri, L, Kn = Ab/At, s, Δs and r

%give each of the variables which change, p,t,a_t->kn,   
    
    dt = mean(diff(t(:,1)));
    Ri = zeros(N,numFires);
    Ro = zeros(N,numFires);
    L  = zeros(N,numFires);
    Ab = zeros(N,numFires); 
    ds = zeros(N,numFires);
    r  = zeros(N,numFires);

    %initial values
    Ri(1,:) = Ri_0;
    Ro(1,:) = Ro_0;
    L(1,:)  = L0;
    % W(1,l)  = Ro(1,l) - Ri(1,l);
    Ab(1,numFires) = Ab_0;
    ds(1,numFires) = 0;
    r(1,numFires)  = 0;
    for l=1:numFires

        Kn(1,l) = Ab_0/A_t(l);
        for k = 1 :(N-1) %solves for ds, regresses grain, updates Ab restart loop
            %our goal is to get regression for each value for pressure to make
            %a curve, Δs -> r so we can then curve fit for a & n
            % fprintf("the iteration is being run")
            % fun = @(ds) ds - (P(k,l)*dt)/(Kn(k,l)*rho*c(l)) %solve for Δs for current step 
            % ds(k,l) = fzero(fun,0); 

            ds(k,l) = (P(k,l)*dt)/(Kn(k,l)*rho*c(l));
            r(k,l) = ds(k,l)/dt; 
            %Regress motor geometry %moooove
            Ri(k+1,l) = Ri(k,l) + ds(k,l);  %inner radius expands outward %r = dels/delt r*dt =ds 
            Ro(k+1,l) = Ro(k,l);  % - r * dt;      % outer radius doesn't regresses inward
            L(k+1,l)  = L(k,l) - 2 * ds(k,l);   % both ends burn inward, turn off if ends inhibited you knew that why am i even 

            Ab(k+1,l)= numGrains*(2*pi*Ri(k+1,l)*L(k+1,l) + 2*pi*(Ro(k+1,l)^2 - Ri(k+1,l)^2));
            Kn(k+1,l) = Ab(k+1,l)/A_t(l);
            W(k+1)  = Ro(k+1,l) - Ri(k+1,l);
            % safeguards
                if Ri(k+1,l) >= Ro(k+1,l) %if Ri would be larger then Ri is set to average between them
                    fprintf("radius safeguard is being run\n") 
                    Ri(k+1,l) = (Ri(k,l) + Ro(k,l)) / (2 - 1e-9) ; 
                end
                if L(k+1,l) < 0
                    fprintf("length safeguard is being run\n")
                    L(k+1,l) = 0;
                end

       end
 end
end

