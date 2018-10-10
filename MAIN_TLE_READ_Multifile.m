%% TWO-LINE-ELEMENT READER / PLOTTER
%
%       Written by:           Tyler Reid (tyreid@alumni.stanford.edu)
%       Lab:                  Stanford GPS Lab
%       Project Start Date:   Oct 09, 2018
%       Last updated:         Oct 09, 2018
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% Given multiple NORAD Two-Line-Element (TLE) files, this matlab code plots
% the orbits of the satellites as well as the Earth. This is meant to be a
% simple orbit plotting / visulization tool. It also plots the Azimuth and
% Elevation angles to these satellites at a given location and time span. 
% 
%% WORKSPACE SETUP

clc
clear 
close all
format longG

% Add path to the Earth plotting function. 
addpath([pwd, '/PlotEarth/PlotEarth']);

% Add path to the TLE data (these are just some examples with a focus on
% GNSS). 
addpath([pwd, '/TLE_Files']);

%% LOAD PHYSICAL CONSTANTS INTO THE GLOBAL WORKSPACE

physical_constants

global mu omega_e R_e

%% USER POSITION & TIME

% Position (lat, long, alt). 
lla = [37.427730, -122.170041, 0];

% Choose a date, by default this chooses the current date/time.
simStart = datenum(clock); % Matlab datenum
simStart = datenum('Sep 29 2018 17:20:00');

% Simulation time. 
tFinal = 6/24; % [days]

%% SELECT THE TLE FILES TO PLOT

% All of GNSS. 
% filenames = {'gps-22'};
% filenames = {'gps-31'};
filenames = {'gps-ops','galileo','beidou','glo-ops'};
% filenames = {'gps-ops'};

% Create colors for plotting. 
colors = lines(length(filenames));

%% DECODE TLE DATA

% Plot the Earth. 
% If you want a color Earth, use 'neompa', 'BlueMarble'.
% If you want a black and white Earth, use 'neomap', 'BlueMarble_bw'.
% A smaller sample step gives a finer resolution Earth.
% h = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 2);

% Compute sidereal time. 
GMST = utc2gmst(datevec(simStart)); % [rad]

% Compute ECEF position of desired location of interest. 
X_POI_ECEF = llh2ECEF(lla(1)*pi/180, lla(2)*pi/180, lla(3));

% Import the TLE data. 
for k = 1:length(filenames)
    % Get the orbital elements. 
    [coe] = two_line_elem_conv(horzcat(filenames{k}, '.txt'), 'all');
    
    % Find latest epoch (all others can be run up from here).
    coeDateNums = datenum(coe.date);
    [val, ind] = min(coeDateNums);
    
    % Define the max time from simStart.
    % a = coe.a(1);
    % n = sqrt(mu/a^3);
    % tFinal = 2*pi/n*1.1/3600/24; % This gives us just over 1 orbit. 
    
    % Create a time vector. 
    tStep = 5; % [sec]
    tSim = simStart:tStep/3600/24:(simStart + tFinal);
    
    % Allocate space.
    RSave = NaN(length(tSim), 3, length(coeDateNums));
    RECEFSave = NaN(length(tSim), 3, length(coeDateNums));
    elSave = NaN(length(tSim), length(coeDateNums));
    azSave = NaN(length(tSim), length(coeDateNums));
    rangeSave = NaN(length(tSim), length(coeDateNums));
    
    % Run through all of the satellites in the TLE file 
    % and compute trajectories for plotting.
    for i = 1:length(coeDateNums) % For each satellite
        for j = 1:length(tSim) % For each time step
            % Get the orbit data. 
            a = coe.a(i);
            n = sqrt(mu / a^3);
            e = coe.e(i);
            inc = coe.i(i) * pi / 180;
            RAAN = coe.RAAN(i) * pi / 180;
            omega = coe.omega(i) * pi / 180;
            M = coe.M(i) * pi / 180 + ...
                n * (tSim(j) - coeDateNums(i)) * 24 * 3600; 
            
            % Convert to ECI and save the data.
            [X_ECI,~] = COE2RV(a, e, inc, RAAN, omega, M);
            RSave(j,:,i) = X_ECI';
            
            % Convert to ECEF. 
            X_ECEF = ECI2ECEF(X_ECI, ...
                GMST + omega_e * (tSim(j) - tSim(1)) * 24 * 3600 );            
            RECEFSave(j,:,i) = X_ECEF; 
            
            % Compute azimuth / elevation. 
            [range,az,el,~,~,~] = rho_az_el(...
                lla(1)*pi/180, lla(2)*pi/180, lla(3), ...
                X_ECEF,[NaN, NaN, NaN]');
            
            % Check if we are below the horizon. 
            if el >= -5 * pi / 180
                rangeSave(j,i) = range; 
                azSave(j,i) = az; 
                elSave(j,i) = el; 
            end
        end
    end
    
    % Plot the orbit (ECI).
    figure(1); 
    for i = 1:length(coeDateNums)
        colorI = k;
        plot3(RSave(:,1,i) / R_e, RSave(:,2,i) / R_e, RSave(:,3,i) / R_e,...
            'color', colors(colorI,:), 'LineWidth', 1)
        plot3(RSave(1,1,i) / R_e, RSave(1,2,i) / R_e, RSave(1,3,i) / R_e,...
            '.', 'color', colors(colorI,:), 'MarkerSize', 10)
        hold on
    end
    
    % Plot the orbit (ECEF). 
    figure(2); 
    for i = 1:length(coeDateNums)
        colorI = k;
        plot3(RECEFSave(:,1,i) / R_e, RECEFSave(:,2,i) / R_e, RECEFSave(:,3,i) / R_e,...
            'color', colors(colorI,:), 'LineWidth', 1)
        plot3(RECEFSave(1,1,i) / R_e, RECEFSave(1,2,i) / R_e, RECEFSave(1,3,i) / R_e,...
            '.', 'color', colors(colorI,:), 'MarkerSize', 10)
        hold on
    end
    
    % Plot elevation as a function of time. 
    figure(3); 
    for i = 1:length(coeDateNums)
        plot(tSim, elSave(:,i) * 180 / pi);
        hold on
    end
    
    % Plot polar plot of satellite visibility. 
    figure(4); 
    for i = 1:length(coeDateNums)
        polar(azSave(:,i), 90 - elSave(:,i) * 180 / pi);
        hold on
    end
    
    % Find start / stop times of RO events and write report to terminal. 
    % TODO make this write to a file. 
    for i = 1:length(coeDateNums)
        % Find zero crossings. 
        x = diff(sign(elSave(:,i)));
        indx_up = find(x>0);
        indx_down = find(x<0);
        
        % Determine which are setting and plot. 
        figure(3); 
        plot(tSim(indx_down), zeros(size(indx_down)), 'o')
        hold on
        
        % Get times, azimuth, and name. 
        if ~isempty(indx_down)
            for j = indx_down
               disp([coe.ID{i},' ', ...
                   coe.PRN{i}, ' ',...
                   coe.designation{i},' ', ...
                   datestr(tSim(j)), ' ', ' Az = ',...
                   num2str(azSave(j,i)*180/pi), ' [deg]' ]); 
            end
        end        
    end
end

%% FORMAT & SAVE FIGURES

figure(1);
axis equal
title('ECI Satellite Trajectories')

figure(2); 
axis equal
title('ECEF Satellite Trajectories')

figure(3); 
ylabel('Elevation [deg]')
xlabel('Date / Time UTC')
grid on
datetick('x')
xlim([simStart, simStart + tFinal])
ylim([-5, 90])
title('Satellite Elevation Angle as a Function of Time')

figure(4); 
view([90 -90])
title('Skyplot')
kill