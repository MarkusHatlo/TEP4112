clear
close all
clc

set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Preambles
GroupName = 'Group 8';
data_filename = '25M';

% Define data paths
data_dir = 'C:\Users\marha\Documents\Skole\TEP4112';
load_dir = [data_dir,'\',GroupName,'\','S_Scans','\'];
save_dir = [data_dir,'\',GroupName,'\processed\'];
precalibration_dir = [data_dir,'\',GroupName,'\Calib\'];
postcalibration_dir = [data_dir,'\',GroupName,'\Calib\'];

savedata = true;
savefilename = [data_filename,'_processed.mat'];

%% Processing
temp_type = 'TP'; % TP for temp. probe; TC for thermocouple
if ~strcmp(temp_type,'TP') && ~strcmp(temp_type,'TC')
    error('Invalid temperature sensor type selection');
end
% Get temperature corrected calibration coefficients
E_gain   =  1;
E_offset = 0;

alfa     = 0.0036;           % Sensor temperature coefficient of resistance
a        = 0.80;             % Overheat ratio
% Constant
R_air   = 286.9;            % Gas constant air
g       = 9.805;            % Gravitational acceleration
C_air   = 1.458e-6;         % Curve fit constant
S_air   = 110.4;            % Curve fit constant

% Calculate gas properties
if exist([precalibration_dir,'CpreHW.mat'],"file") && exist([postcalibration_dir,'CpostHW.mat'],"file")
    calib_type = 3;
elseif exist([precalibration_dir,'CpostHW.mat'],"file")
    calib_type = 2;
elseif exist([postcalibration_dir,'CpreHW.mat'],"file")
    calib_type = 1;
else
    error('No calibration file exists');
end

if calib_type == 3
    %%% Pre Calib %%%
    load([precalibration_dir,'CpreHW.mat']);
    timestamp_pre = (datenum(CHW.start_time)+datenum(CHW.end_time))/2;
    if strcmp(temp_type,'TP')
        T = CHW.T_tp;
    else
        T = CHW.T_tc;
    end
    Patm = CHW.P_atm;
    P_atm = CHW.P_atm;                % Keeping end pressure for hotwire process
    sgHg    = 13.6 - 0.0024*T;          % Spesific gravity Hg
    p_atm   = sgHg*Patm*g;                  % Atmospheric pressure in Pa
    % Calculate wire temperature
    T_wire  = CHW.T_ref + a/alfa;
    E = CHW.SW_E/E_gain + E_offset;
    % Calculate density and viscosity
    T_K  = T + 273.15;                        % Temperature in K
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    k       = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));
    
    % Calculate temperature difference
    delta_T = T_wire - T_K;
    % Find curve fit constants
    C_pre = polyfit((E.^2)./(k.*delta_T),CHW.U_inf./nu_f,4);
    
    % Generate points to plot polyfit function
    minX = min((E.^2)./(k.*delta_T));
    maxX = max((E.^2)./(k.*delta_T));
    X_pol   = linspace(minX,maxX);
    Y_pol   = C_pre(1)*X_pol.^4 + C_pre(2)*X_pol.^3 + C_pre(3)*X_pol.^2 + C_pre(4)*X_pol + C_pre(5);
    
    % Plot polyfit curves
    figure()
    hold on
    plot(X_pol,Y_pol,'b',(E.^2)./(k.*delta_T),CHW.U_inf./nu_f,'bo')
    grid on
    xlabel('$E^2/(k*\Delta T)$'); ylabel('$U/\nu$')
    
    %%% Post Calib %%%
    load([postcalibration_dir,'CpostHW.mat']);
    timestamp_post = (datenum(CHW.start_time)+datenum(CHW.end_time))/2;
    if strcmp(temp_type,'TP')
        T = CHW.T_tp;
    else
        T = CHW.T_tc;
    end
    Patm = CHW.P_atm;
    P_atm = (P_atm + CHW.P_atm)/2;      % Averaging pre-end- and post-start-pressure
    sgHg    = 13.6 - 0.0024*T;          % Spesific gravity Hg
    p_atm   = sgHg*Patm*g;                  % Atmospheric pressure in Pa
    % Calculate wire temperature
    T_wire  = CHW.T_ref + a/alfa;
    E = CHW.SW_E/E_gain + E_offset;
    % Calculate density and viscosity
    T_K  = T + 273.15;                        % Temperature in K
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    k       = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));
    
    % Calculate temperature difference
    delta_T = T_wire - T_K;
    % Find curve fit constants
    C_post = polyfit((E.^2)./(k.*delta_T),CHW.U_inf./nu_f,4);
    
        % Generate points to plot polyfit function
        minX = min((E.^2)./(k.*delta_T));
        maxX = max((E.^2)./(k.*delta_T));
        X_pol   = linspace(minX,maxX);
        Y_pol   = C_post(1)*X_pol.^4 + C_post(2)*X_pol.^3 + C_post(3)*X_pol.^2 + C_post(4)*X_pol + C_post(5);
    
        % Plot polyfit curves
        plot(X_pol,Y_pol,'r',(E.^2)./(k.*delta_T),CHW.U_inf./nu_f,'ro')
        legend('Pre-Calib fit','Pre-Calib pts','Post-Calib fit','Post-Calib pts','Location','northwest')
    
    %%% Load HW data %%%
    load([load_dir,data_filename,'.mat']);
    if strcmp(temp_type,'TP')
        T_K = HWA.T_tp_RAW + 273.15;
    else
        try
            T_K = linspace(mean(HWA.T_tc_start_RAW),mean(HWA.T_tc_end_RAW),numel(HWA.SW_RAW))' + 273.15;
        catch
            T_K = mean(HWA.raw_T_tc) + 273.15;
        end
    end
    
    sgHg    = 13.6 - 0.0024*(T_K-273.15);          % Spesific gravity Hg
    try
        p_atm   = sgHg*P_atm*g;                  % Atmospheric pressure in Pa
    catch
        p_atm   = sgHg*Patm*g;                  % Atmospheric pressure in Pa
    end
    
    % Gas properties
    rho     = p_atm./(R_air.*T_K);             % Density of the air [kg/m^3]
    mu      = C_air*T_K.^(3/2)./(S_air+T_K);  % Dynamic viscosity
    nu      = mu./rho;
    
    T_wire = CHW.T_ref + a/alfa;
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    
    E_HW    = HWA.SW_RAW/E_gain + E_offset;
    k = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));
    cal_for  = (E_HW.^2)./(k.*(T_wire-T_K));
    u_pre = nu_f.*(C_pre(1)*cal_for.^4 + C_pre(2)*cal_for.^3 + C_pre(3)*cal_for.^2 + C_pre(4)*cal_for + C_pre(5));
    u_post = nu_f.*(C_post(1)*cal_for.^4 + C_post(2)*cal_for.^3 + C_post(3)*cal_for.^2 + C_post(4)*cal_for + C_post(5));
    timestamp = datenum((HWA.time_start+HWA.time_end)/2);
    u = u_pre + (timestamp - timestamp_pre)*(u_post-u_pre)/(timestamp_post-timestamp_pre);
    
elseif calib_type == 2
    %%% Post Calib %%%
    load([postcalibration_dir,'CpostHW.mat']);
    if strcmp(temp_type,'TP')
        T = CHW.T_tp;
    else
        T = CHW.T_tc;
    end
    Patm = (CHW.P_atm_1 + CHW.P_atm_2)/2;
    sgHg    = 13.6 - 0.0024*T;          % Spesific gravity Hg
    p_atm   = sgHg*Patm*g;                  % Atmospheric pressure in Pa
    
    % Calculate density and viscosity
    T_K  = T + 273.15;                        % Temperature in K
        % Calculate wire temperature
    T_wire  = CHW.T_ref + a/alfa;
    E = CHW.SW_E/E_gain + E_offset;
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    k       = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));
    
    % Calculate temperature difference
    delta_T = T_wire - T_K;
    % Find curve fit constants
    C = polyfit((E.^2)./(k.*delta_T),CHW.U_inf./nu_f,4);
    
    %%% Load HW data %%%
    load([load_dir,data_filename,'.mat']);
    if strcmp(temp_type,'TP')
        T_K = HWA.T_tp_RAW + 273.15;
    else
        T_K = linspace(mean(HWA.T_tc_start_RAW),mean(HWA.T_tc_end_RAW),numel(HWA.SW_RAW))' + 273.15;
    end
    
    sgHg    = 13.6 - 0.0024*(T_K-273.15);          % Spesific gravity Hg
    p_atm   = sgHg*HWA.P_atm*g;                  % Atmospheric pressure in Pa
    
    % Gas properties
    rho     = p_atm./(R_air.*T_K);             % Density of the air [kg/m^3]
    mu      = C_air*T_K.^(3/2)./(S_air+T_K);  % Dynamic viscosity
    nu      = mu./rho;
    
    T_wire = CHW.T_ref + a/alfa;
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    
    E_HW    = HWA.SW_RAW/E_gain + E_offset;
    k = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));
    cal_for  = (E_HW.^2)./(k.*(T_wire-T_K));
    u = nu_f.*(C(1)*cal_for.^4 + C(2)*cal_for.^3 + C(3)*cal_for.^2 + C(4)*cal_for + C(5));

elseif calib_type == 1
    %%% Pre Calib %%%
    load([precalibration_dir,'CpreHW.mat']);
    if strcmp(temp_type,'TP')
        T = CHW.T_tp;
    else
        T = CHW.T_tc;
    end
    Patm = (CHW.P_atm_1 + CHW.P_atm_2)/2;
    sgHg    = 13.6 - 0.0024*T;          % Spesific gravity Hg
    p_atm   = sgHg*Patm*g;                  % Atmospheric pressure in Pa
    % Calculate wire temperature
    T_wire  = CHW.T_ref + a/alfa;
    E = CHW.SW_E/E_gain + E_offset;
    % Calculate density and viscosity
    T_K  = T + 273.15;                        % Temperature in K
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    k       = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));    
    
    % Calculate temperature difference
    delta_T = T_wire - T_K;
    % Find curve fit constants
    C = polyfit((E.^2)./(k.*delta_T),CHW.U_inf./nu_f,4);
    
    %%% Load HW data %%%
    load([load_dir,data_filename,'.mat']);
    if strcmp(temp_type,'TP')
        T_K = HWA.T_tp_RAW + 273.15;
    else
        T_K = linspace(mean(HWA.T_tc_start_RAW),mean(HWA.T_tc_end_RAW),numel(HWA.SW_RAW))' + 273.15;
    end
    
    sgHg    = 13.6 - 0.0024*(T_K-273.15);          % Spesific gravity Hg
    p_atm   = sgHg*HWA.P_atm*g;                  % Atmospheric pressure in Pa
    
    % Gas properties
    rho     = p_atm./(R_air.*T_K);             % Density of the air [kg/m^3]
    mu      = C_air*T_K.^(3/2)./(S_air+T_K);  % Dynamic viscosity
    nu      = mu./rho;
    
    T_wire = CHW.T_ref + a/alfa;
    Tf_K = (T_wire+T_K)/2; % air film temperature
    rho_f     = p_atm./(R_air.*Tf_K);               % Density of the air
    mu_f      = C_air*Tf_K.^(3/2)./(S_air+Tf_K);  % Dynamic viscosity
    nu_f      = mu_f./rho_f;                              % Kinematic viscosity
    
    E_HW    = HWA.SW_RAW/E_gain + E_offset;
    k = 418.4*(5.75*10^-5*(1+0.00317*Tf_K - 0.0000021*Tf_K.^2));
    cal_for  = (E_HW.^2)./(k.*(T_wire-T_K));
    u = nu_f.*(C(1)*cal_for.^4 + C(2)*cal_for.^3 + C(3)*cal_for.^2 + C(4)*cal_for + C(5));    
end

rho = mean(rho);
T = mean(T_K) - 273.15;
nu = mean(nu);
fs = HWA.fs_actual;
ts = HWA.ts_actual;

%% Plotting the spectrum

FONT_SIZE  = 16;
FONT_SIZE2 = 10;

figure()
plot([1:numel(u)]'/fs,u);
set(gca,'fontsize',FONT_SIZE2);
xl = xlabel('$t$ [s]');
set(xl,'unit','character');
set(xl,'interpreter','latex');
set(xl,'fontsize', FONT_SIZE);
yl = ylabel('$u_1$ [m/s]');
set(yl,'unit','character');
set(yl,'interpreter','latex');
set(yl,'fontsize', FONT_SIZE);
grid on
box on

%% Saving data
if savedata
    disp('Saving processed data...');
    % Create the output folder if it doesn´t exist
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
    % save([save_dir savefilename]);
    save([save_dir savefilename],'u','P_atm','T','rho','nu','fs','ts','-v6');
end