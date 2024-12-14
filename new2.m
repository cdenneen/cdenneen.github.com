% Constants
solar_constant = 1353; % W/m²
Cn = 0.75; % PV system efficiency
latitude = 42.984333; % Syracuse, NY
longitude = -76.142167;
tilt_angle = 30; % Collector tilt angle (degrees)
surface_azimuth = 10; % Collector azimuth angle (degrees)
albedo = 0.2; % Ground reflectance
household_consumption = 10000; % Yearly household electricity consumption in kWh

% Dates for simulation (one day per month)
dates = ["2024-01-15", "2024-02-15", "2024-03-15", "2024-04-15", "2024-05-15", "2024-06-15", "2024-07-15", "2024-08-15", "2024-09-15", "2024-10-15", "2024-11-15", "2024-12-15"];

% Times from 4 AM to 9 PM (in 1-hour increments)
times = arrayfun(@(x) sprintf('%02d:00', x), 4:21, 'UniformOutput', false);

% Calculate daily solar radiation for each date
daily_radiation = zeros(1, length(dates));
for i = 1:length(dates)
    date = dates{i};
    radiation_values = zeros(1, length(times));
    for j = 1:length(times)
        time = times{j};
        Ic = solar_radiation(latitude, longitude, date, time, tilt_angle, surface_azimuth, albedo, Cn, 0.144, 0.06);
        radiation_values(j) = Ic;
    end
    daily_radiation(i) = sum(radiation_values);
end

% Calculate yearly radiation by averaging daily values and scaling
average_daily_radiation = mean(daily_radiation); % W/m²
total_yearly_radiation = average_daily_radiation * 365; % W/m² per year

% Convert to kWh/m² per year
total_yearly_radiation_kWh = total_yearly_radiation / 1000; % kWh/m² per year

% Calculate the required PV area
required_PV_area = household_consumption / total_yearly_radiation_kWh; % m²

fprintf('Required PV area to meet yearly household electricity demand: %.2f m²\n', required_PV_area);

% Define the function to calculate solar radiation
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