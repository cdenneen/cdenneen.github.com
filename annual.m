% Constants
solar_constant = 1353; % W/m²
Cn = 0.75; % PV system efficiency
latitude = 42.984333; % Syracuse, NY
longitude = -76.142167;
tilt_angle = 30; % Assumed tilt angle (degrees)
surface_azimuth = 10; % Assumed collector azimuth angle (degrees)
albedo = 0.2; % Ground reflectance
household_consumption = 10000; % Yearly household electricity consumption in kWh
conversion_efficiency = 0.15; % 15% solar to electricity conversion
electricity_cost_per_kWh = 0.13; % Cost of electricity per kWh (USD)

% Dates for simulation (one day per month)
dates = ["2024-01-15", "2024-02-15", "2024-03-15", "2024-04-15", "2024-05-15", "2024-06-15", "2024-07-15", "2024-08-15", "2024-09-15", "2024-10-15", "2024-11-15", "2024-12-15"];
labels = ["January 15", "February 15", "March 15", "April 15", "May 15", "June 15", "July 15", "August 15", "September 15", "October 15", "November 15", "December 15"];

% Times from 4 AM to 9 PM (in 1-hour increments)
times = arrayfun(@(x) sprintf('%02d:00', x), 4:21, 'UniformOutput', false);

% Variables to store monthly energy and total yearly energy
monthly_energy = zeros(1, length(dates));
total_yearly_energy = 0;

for i = 1:length(dates)
    date = dates{i};
    radiation_values = zeros(1, length(times));

    for j = 1:length(times)
        time = times{j};
        Ic = solar_radiation(latitude, longitude, date, time, tilt_angle, surface_azimuth, albedo, Cn, 0.144, 0.06, solar_constant);
        radiation_values(j) = Ic;
    end

    % Calculate daily energy in kWh/m²
    daily_energy = sum(radiation_values) / 1000; % Convert Wh/m² to kWh/m²
    monthly_energy(i) = daily_energy * 30; % Approximate monthly energy (30 days)
end

% Sum up total annual energy
total_yearly_energy = sum(monthly_energy);

% Calculate electricity generated by the panel
electricity_generated = total_yearly_energy * conversion_efficiency; % kWh/m²

% Calculate the required PV area
required_PV_area = household_consumption / electricity_generated; % m²

% Calculate the annual electricity savings cost
annual_savings_cost = electricity_generated * required_PV_area * electricity_cost_per_kWh; % USD

% Display results
fprintf('\nMonthly Energy (kWh/m²):\n');
for i = 1:length(dates)
    fprintf('%s: %.2f kWh/m²\n', labels{i}, monthly_energy(i));
end
fprintf('\nTotal Annual Energy: %.2f kWh/m²\n', total_yearly_energy);
fprintf('Electricity Generated by the Panel: %.2f kWh/m²\n', electricity_generated);
fprintf('Required PV Area to meet household consumption: %.2f m²\n', required_PV_area);
fprintf('Annual Electricity Savings Cost: $%.2f\n', annual_savings_cost);

% Define the function to calculate solar radiation
function Ic = solar_radiation(lat, lon, date, time, tilt_angle, surface_azimuth, albedo, Cn, k, C, solar_constant)
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