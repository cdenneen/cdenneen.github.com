% Solar Radiation Calculation and Plotting
% Created on Mon Nov 18 23:21:19 2024
% Author: Kalhaku Mclester

% Constants
solar_constant = 1353; % W/m²

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
                (cos(declination) * cos(lat) * cos(tilt_angle) * cos(solar_hour_angle))