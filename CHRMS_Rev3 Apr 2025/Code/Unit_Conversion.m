%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit_Conversion
% 2025/03/23
% Adrian Comisso
%
% Description: 
% This function converts values between metric and non-metric units
% based on the unit category (pressure, temperature, density, etc.)
% The function infers the unit category from the unit string and applies
% the appropriate conversion.
%
% Inputs:
% value - the numerical value to convert
% unit - string representing the unit (e.g., 'kg/m^2', 'K', 'Pa', etc.)
% to_metric - boolean (true: convert to metric, false: convert from metric)
%
% Outputs:
% result - the converted value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = Unit_Conversion(value, unit, to_metric)
    % Determine the unit category
    category = determineCategory(unit);
    
    % Apply the appropriate conversion based on the category
    switch category
        case 'pressure'
            result = convertPressure(value, unit, to_metric);
        case 'temperature'
            result = convertTemperature(value, unit, to_metric);
        case 'density'
            result = convertDensity(value, unit, to_metric);
        case 'volume'
            result = convertVolume(value, unit, to_metric);
        case 'mass'
            result = convertMass(value, unit, to_metric);
        case 'length'
            result = convertLength(value, unit, to_metric);
        case 'percent'
            result = convertPercent(value, unit, to_metric);
        otherwise
            error('Unknown unit category for unit: %s', unit);
    end
end

% Helper function to determine the category of a unit
function [category] = determineCategory(unit)
    % Pressure units
    if any(strcmp(unit, {'Pa', 'kPa', 'MPa', 'Atm', 'Bar', 'psi'}))
        category = 'pressure';
    % Temperature units
    elseif any(strcmp(unit, {'K', 'C', 'R', 'F'}))
        category = 'temperature';
    % Density units
    elseif any(strcmp(unit, {'kg/m^3', 'g/cm^3', 'lbm/in^3', 'lbm/ft^3'}))
        category = 'density';
    % Volume units
    elseif any(strcmp(unit, {'m^3', 'cm^3', 'L', 'Gal', 'in^3', 'ft^3'}))
        category = 'volume';
    % Mass units
    elseif any(strcmp(unit, {'kg', 'g', 'oz', 'lbm', 'psf'}))
        category = 'mass';
    % Length units
    elseif any(strcmp(unit, {'m', 'cm', 'mm', 'in', 'ft'}))
        category = 'length';
    % Percent units
    elseif any(strcmp(unit, {'%', 'dec'}))
        category = 'percent';
    else
        error('Unknown unit: %s', unit);
    end
end

% Pressure conversion function
function [result] = convertPressure(value, unit, to_metric)
    % Base metric unit for pressure is Pa
    if to_metric
        % Convert from given unit to Pa
        switch unit
            case 'Pa'
                result = value; % Already in Pa
            case 'kPa'
                result = value * 1000; % kPa to Pa
            case 'MPa'
                result = value * 1e6; % MPa to Pa
            case 'Atm'
                result = value * 101325; % Atm to Pa
            case 'Bar'
                result = value * 1e5; % Bar to Pa
            case 'psi'
                result = value * 6894.76; % psi to Pa
            otherwise
                error('Unknown pressure unit: %s', unit);
        end
    else
        % Convert from Pa to given unit
        switch unit
            case 'Pa'
                result = value; % Already in Pa
            case 'kPa'
                result = value / 1000; % Pa to kPa
            case 'MPa'
                result = value / 1e6; % Pa to MPa
            case 'Atm'
                result = value / 101325; % Pa to Atm
            case 'Bar'
                result = value / 1e5; % Pa to Bar
            case 'psi'
                result = value / 6894.76; % Pa to psi
            otherwise
                error('Unknown pressure unit: %s', unit);
        end
    end
end

% Temperature conversion function
function [result] = convertTemperature(value, unit, to_metric)
    % Base metric unit for temperature is K
    if to_metric
        % Convert from given unit to K
        switch unit
            case 'K'
                result = value; % Already in K
            case 'C'
                result = value + 273.15; % C to K
            case 'R'
                result = value * 5/9; % R to K
            case 'F'
                result = (value + 459.67) * 5/9; % F to K
            otherwise
                error('Unknown temperature unit: %s', unit);
        end
    else
        % Convert from K to given unit
        switch unit
            case 'K'
                result = value; % Already in K
            case 'C'
                result = value - 273.15; % K to C
            case 'R'
                result = value * 9/5; % K to R
            case 'F'
                result = value * 9/5 - 459.67; % K to F
            otherwise
                error('Unknown temperature unit: %s', unit);
        end
    end
end

% Density conversion function
function [result] = convertDensity(value, unit, to_metric)
    % Base metric unit for density is kg/m^3
    if to_metric
        % Convert from given unit to kg/m^3
        switch unit
            case 'kg/m^3'
                result = value; % Already in kg/m^3
            case 'g/cm^3'
                result = value * 1000; % g/cm^3 to kg/m^3
            case 'lbm/in^3'
                result = value * 27679.9; % lbm/in^3 to kg/m^3
            case 'lbm/ft^3'
                result = value * 16.0185; % lbm/ft^3 to kg/m^3
            otherwise
                error('Unknown density unit: %s', unit);
        end
    else
        % Convert from kg/m^3 to given unit
        switch unit
            case 'kg/m^3'
                result = value; % Already in kg/m^3
            case 'g/cm^3'
                result = value / 1000; % kg/m^3 to g/cm^3
            case 'lbm/in^3'
                result = value / 27679.9; % kg/m^3 to lbm/in^3
            case 'lbm/ft^3'
                result = value / 16.0185; % kg/m^3 to lbm/ft^3
            otherwise
                error('Unknown density unit: %s', unit);
        end
    end
end

% Volume conversion function
function [result] = convertVolume(value, unit, to_metric)
    % Base metric unit for volume is m^3
    if to_metric
        % Convert from given unit to m^3
        switch unit
            case 'm^3'
                result = value; % Already in m^3
            case 'cm^3'
                result = value / 1e6; % cm^3 to m^3
            case 'L'
                result = value / 1000; % L to m^3
            case 'Gal'
                result = value * 0.00378541; % Gal to m^3
            case 'in^3'
                result = value * 1.63871e-5; % in^3 to m^3
            case 'ft^3'
                result = value * 0.0283168; % ft^3 to m^3
            otherwise
                error('Unknown volume unit: %s', unit);
        end
    else
        % Convert from m^3 to given unit
        switch unit
            case 'm^3'
                result = value; % Already in m^3
            case 'cm^3'
                result = value * 1e6; % m^3 to cm^3
            case 'L'
                result = value * 1000; % m^3 to L
            case 'Gal'
                result = value / 0.00378541; % m^3 to Gal
            case 'in^3'
                result = value / 1.63871e-5; % m^3 to in^3
            case 'ft^3'
                result = value / 0.0283168; % m^3 to ft^3
            otherwise
                error('Unknown volume unit: %s', unit);
        end
    end
end

% Mass conversion function
function [result] = convertMass(value, unit, to_metric)
    % Base metric unit for mass is kg
    if to_metric
        % Convert from given unit to kg
        switch unit
            case 'kg'
                result = value; % Already in kg
            case 'g'
                result = value / 1000; % g to kg
            case 'oz'
                result = value * 0.0283495; % oz to kg
            case 'lbm'
                result = value * 0.453592; % lbm to kg
            case 'psf'
                result = value * 4.88243; % psf to kg (assuming pound-force per square foot)
            otherwise
                error('Unknown mass unit: %s', unit);
        end
    else
        % Convert from kg to given unit
        switch unit
            case 'kg'
                result = value; % Already in kg
            case 'g'
                result = value * 1000; % kg to g
            case 'oz'
                result = value / 0.0283495; % kg to oz
            case 'lbm'
                result = value / 0.453592; % kg to lbm
            case 'psf'
                result = value / 4.88243; % kg to psf (assuming pound-force per square foot)
            otherwise
                error('Unknown mass unit: %s', unit);
        end
    end
end

% Length conversion function
function [result] = convertLength(value, unit, to_metric)
    % Base metric unit for length is m
    if to_metric
        % Convert from given unit to m
        switch unit
            case 'm'
                result = value; % Already in m
            case 'cm'
                result = value / 100; % cm to m
            case 'mm'
                result = value / 1000; % mm to m
            case 'in'
                result = value * 0.0254; % in to m
            case 'ft'
                result = value * 0.3048; % ft to m
            otherwise
                error('Unknown length unit: %s', unit);
        end
    else
        % Convert from m to given unit
        switch unit
            case 'm'
                result = value; % Already in m
            case 'cm'
                result = value * 100; % m to cm
            case 'mm'
                result = value * 1000; % m to mm
            case 'in'
                result = value / 0.0254; % m to in
            case 'ft'
                result = value / 0.3048; % m to ft
            otherwise
                error('Unknown length unit: %s', unit);
        end
    end
end

% Percent conversion function
function [result] = convertPercent(value, unit, to_metric)
    % Base metric unit for percent is decimal (dec)
    if to_metric
        % Convert from given unit to decimal
        switch unit
            case 'dec'
                result = value; % Already in decimal
            case '%'
                result = value / 100; % % to decimal
            otherwise
                error('Unknown percent unit: %s', unit);
        end
    else
        % Convert from decimal to given unit
        switch unit
            case 'dec'
                result = value; % Already in decimal
            case '%'
                result = value * 100; % decimal to %
            otherwise
                error('Unknown percent unit: %s', unit);
        end
    end
end