% Solar Radiation Calculation and Plotting
% Created on Mon Nov 18 23:21:19 2024
% Author: Kalhaku Mclester

% Constants
solar_constant = 1353; % W/m²

% Parameters
latitude = 42.984333; % Syracuse, NY
longitude = -76.142167;
tilt_angle = 30; % Collector tilt angle (degrees)
surface_azimuth = 10; % Collector azimuth angle (degrees)
albedo = 0.2; % Ground reflectance

% Dates for different seasons
dates = ["2024-03-21", "2024-06-21", "2024-09-21", "2024-12-21"];
labels = ["March 21", "June 21", "September 21", "December 21"];

% Times from 4 AM to 9 PM (in 1-minute increments)
times = arrayfun(@(x) sprintf('%02d:%02d', floor(x/60), mod(x, 60)), 4*60:21*60-1, 'UniformOutput', false);

% Sunrise and Sunset Errors (in HH:MM format)
actual_times = containers.Map({'March 21', 'June 21', 'September 21', 'December 21'}, ...
    {{'7:05 AM', '7:18 PM'}, {'5:25 AM', '8:47 PM'}, {'6:51 AM', '7:03 PM'}, {'7:33 AM', '4:33 PM'}});

% Plotting
figure('Position', [100, 100, 1200, 700]);
hold on;

% Calculate and plot solar radiation
for i = 1:length(dates)
    date = dates{i};
    label = labels{i};
    radiation_values = zeros(1, length(times));
    non_zero_times = {};

    for j = 1:length(times)
        time = times{j};
        Ic = solar_radiation(latitude, longitude, date, time, tilt_angle, surface_azimuth, albedo, 1, 0.144, 0.06);
        radiation_values(j) = Ic;
        if Ic > 0
            non_zero_times{end+1} = time;
        end
    end

    % Plot solar radiation
    plot(1:length(times), radiation_values, 'DisplayName', label);

    % Calculate sunrise and sunset
    if ~isempty(non_zero_times)
        first_non_zero = datetime(non_zero_times{1}, 'InputFormat', 'HH:mm') + hours(1);
        last_non_zero = datetime(non_zero_times{end}, 'InputFormat', 'HH:mm') + hours(1);
        predicted_sunrise = datestr(first_non_zero, 'HH:MM PM');
        predicted_sunset = datestr(last_non_zero, 'HH:MM PM');
    else
        predicted_sunrise = 'N/A';
        predicted_sunset = 'N/A';
    end

    % Print the times and errors
    actual_sunrise = datetime(actual_times(label){1}, 'InputFormat', 'hh:mm a');
    actual_sunset = datetime(actual_times(label){2}, 'InputFormat', 'hh:mm a');

    if ~isempty(non_zero_times)
        sunrise_error = minutes(first_non_zero - actual_sunrise);
        sunset_error = minutes(last_non_zero - actual_sunset);
    else
        sunrise_error = NaN;
        sunset_error = NaN;
    end

    fprintf('%s:\n', label);
    fprintf('  Predicted Sunrise: %s, Actual Sunrise: %s, Error: %.2f minutes\n', predicted_sunrise, actual_times(label){1}, sunrise_error);
    fprintf('  Predicted Sunset: %s, Actual Sunset: %s, Error: %.2f minutes\n', predicted_sunset, actual_times(label){2}, sunset_error);
end

% Plot settings
xlabel('Time');
ylabel('Solar Radiation (W/m²)');
title('Solar Radiation throughout the Day for Different Seasons');
legend('show');
grid on;
hold off;

% Function to calculate solar radiation
function Ic = solar_radiation(lat, lon, date, time, tilt_angle, surface_azimuth, albedo, Cn, k, C)
    % Convert inputs to radians
    lat = deg2rad(lat);
    tilt_angle = deg2rad(tilt_angle);
    surface_azimuth = deg2rad(surface_azimuth);

    % Parse the date and time
    datetime_obj = datetime(date + " " + time, 'InputFormat', 'yyyy-MM-dd HH:mm');

    % Calculate day of the year (n)
    day_of_year = day(datetime_obj, 'dayofyear');

    % Equation of Time (ET)
    B = deg2rad(360 * (day_of_year - 81) / 364); % Radians
    ET = 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B); % Minutes

    % Convert EST to Solar Time
    local_standard_meridian = -75; % For EST (degrees)
    solar_time_offset = ET + 4 * (lon - local_standard_meridian);
    solar_time = datetime_obj + minutes(solar_time_offset);
    solar_hour_angle = deg2rad(15 * (hour(solar_time) + minute(solar_time) / 60 - 12)); % Radians

    % Solar Declination (delta)
    declination = deg2rad(23.45 * sin(deg2rad(360 / 365 * (day_of_year - 81))));

    % Solar Altitude Angle (alpha)
    sin_alpha = sin(lat) * sin(declination) + cos(lat) * cos(declination) * cos(solar_hour_angle);
    alpha = asin(sin_alpha);

    % Direct Normal Irradiance (Ib,n)
    if rad2deg(alpha) > 0
        Ib_n = Cn * solar_constant * (1 + 0.034 * cos(deg2rad(360 * day_of_year / 365.25))) * exp(-k / sin(alpha));
    else
        Ib_n = 0;
    end

    % Diffuse Radiation (Id)
    Id = C * Ib_n * cos(tilt_angle / 2) * cos(tilt_angle / 2);

    % Ground Reflected Radiation (Ir)
    Ir = albedo * Ib_n * sin(alpha + C) * sin(tilt_angle / 2) * sin(tilt_angle / 2);

    % Angle of Incidence (theta)
    cos_theta = (sin(declination) * sin(lat) * cos(tilt_angle)) - ...
                (sin(declination) * cos(lat) * sin(tilt_angle) * cos(surface_azimuth)) + ...
                (cos(declination) * cos(lat) * cos(tilt_angle) * cos(solar_hour_angle)) + ...
                (cos(declination) * sin(lat) * sin(tilt_angle) * cos(surface_azimuth) * cos(solar_hour_angle)) + ...
                (cos(declination) * sin(tilt_angle) * sin(surface_azimuth) * sin(solar_hour_angle));
    cos_theta = max(cos_theta, 0);

    % Beam Radiation on the Collector (Ib_n * cos(theta))
    Ib_cos_theta = Ib_n * cos_theta;

    % Total Radiation on the Collector (Ic)
    Ic = Ib_cos_theta + Id + Ir;
end