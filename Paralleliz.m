%% define constant
% permeability
meu = 12.56637061e-7;
% permittivity
epsln_0 = 8.854e-12;
% signal frequency
f = 970e6;
%signal speed
c = (1/sqrt(meu * epsln_0));
omega = 2 * pi * f;
% eta and beta
bt = abs(2 * pi * f * sqrt(meu * epsln_0));  
et = sqrt(meu / epsln_0);  

% qw: Quarter wavelength

qw = (c / f) / 4; 
hw=(c / f) / 2; 
%% source data
s_x = 0;
s_y = 442.0;
s_ht = 442.0 - 390.0;

%% load the terrian date from file
fileID = fopen("X.04", 'r');
formatSpec = "%d %f";
sizeA = [2 inf];
A = fscanf(fileID, formatSpec,sizeA);
A = A';
% terrian space coordiation
x = A(:,1);
y = A(:,2);
% group size
grpsz = 13;
% distance between two points 
grpstp = 10.0;
% terrain length
terr_length = 700;
% Total group
tot_groups = floor(terr_length /(grpsz * qw));
tot_groups_hw = floor(terr_length /(grpsz * hw));

surface_current = zeros(tot_groups+1, 1);
surface_current_hw = zeros(tot_groups_hw+1, 1);

electric_field = zeros(tot_groups+1, 1);
electric_field_hw = zeros(tot_groups_hw+1, 1);

ef_db = zeros(tot_groups+1, 2);
ef_db_hw = zeros(tot_groups_hw+1, 2);

surface_abscurrent = zeros(tot_groups+1, 2);
surface_abscurrent_hw = zeros(tot_groups_hw+1, 2);

sum = complex(0.0, 0.0);
sum_hw = complex(0.0, 0.0);

for n = 1:grpsz-1
    sum = sum + abs(H02((grpsz-n) * qw * bt));
    sum_hw = sum + abs(H02((grpsz-n) * hw * bt));
end
const = sum + abs(complex(1.0, 1.05)); 
const_hw = sum_hw + abs(complex(1.0, 1.05)); 
zself= complex(1.0, 1.05);
zself_hw = complex(1.0, 1.05);
% Surface current value
for n = 1:tot_groups+1
    sum_current = complex(0.0, 0.0);
    
    for m = 95:(n-1)
        sum_current = sum_current + const * surface_current(m) * H02(bt * R_p_q(m, n, qw, y));
    end
    surface_current(n) = (Ei_Rad(bt *R_source_p(n, qw, y))- sum_current) /zself;
end

for n_hw = 1:tot_groups_hw+1
    sum_current_hw = complex(0.0, 0.0);
    for m_hw = 95:(n_hw-1)
     sum_current_hw = sum_current_hw + const_hw * surface_current(m_hw) * H02(bt * R_p_q_hw(m_hw, n_hw, hw, y));
        
    end
    surface_current_hw(n_hw) = (Ei_Rad(bt * R_source_p_hw(n_hw, hw, y)) - sum_current_hw) /zself;
end


for n = 1:tot_groups+1
    electric_field(n) = Ei_Rad(bt * R_source_obs(n, qw, y));
    
    for m = 95:(n-1)
        electric_field(n) = electric_field(n) - const * surface_current(m) * H02(bt * R_surface_obs(m, n, qw, y));
        
    end
    ef_db(n,1) = X(n, qw);
    
    ef_db(n,2) = 10 * log10(abs(electric_field(n))/sqrt(R_source_obs(n+1, qw, y)));
   
end


for k = 1:tot_groups_hw+1
   
    electric_field_hw(k) = Ei_Rad(bt * R_source_obs_hw(k, hw, y));
    for q = 95:(n_hw-1)
        
        electric_field_hw(k) = electric_field_hw(k) - const_hw * surface_current_hw(q) * H02(bt * R_surface_obs_hw(q, k, hw, y));
    end
    
     ef_db_hw(k,1) = X(k, hw);
    
    ef_db_hw(k,2) = 10 * log10(abs(electric_field_hw(k))/sqrt(R_source_obs_hw(k+1, hw, y)));
end



%% plots
f1 = figure;
plot(ef_db(:,1), ef_db(:,2));
title("Electric Field for Quarter Wavelength Reception");
ylabel("Electric Field (E) in dB");
xlabel("Distance(X)in Meters");

f2= figure;
plot(ef_db_hw(:,1), smoothdata(ef_db_hw(:,2)));
title("Electric Field for Half Wavelength Reception");
ylabel("Electric Field (E) in dB");
xlabel("Distance(X)in Meters");


for n = 1:tot_groups+1
    surface_abscurrent(n, 1) = X(n, qw);    
    surface_abscurrent(n, 2) = abs(surface_current(n));
end


for n_hw = 1:tot_groups_hw+1
   surface_abscurrent_hw(n_hw, 1) = X(n_hw, hw);    
     surface_abscurrent_hw(n_hw, 2) = abs(surface_current_hw(n_hw));
end


f3 = figure;
plot(surface_abscurrent(:,1), surface_abscurrent(:,2));
title("Surface Current Value");
ylabel("Surface Current(J)");
xlabel("Distance(X)in Meters");


f4 = figure;
plot(surface_abscurrent_hw(:,1),smoothdata(surface_abscurrent_hw(:,2)));
title("Surface Current Value");
ylabel("Surface Current(J)");
xlabel("Distance(X)in Meters");

% Calculate power
power=zeros(tot_groups+1,1);
power_hw=zeros(tot_groups_hw+1,1);

for z = 1:tot_groups+1
       power(z,1) = surface_abscurrent(z,2)* abs(electric_field(z,1));  
end

for z_hw = 1:tot_groups_hw+1        
    power_hw(z_hw,1) = surface_abscurrent_hw(z_hw,2)* abs(electric_field_hw(z_hw,1)); 
end

powerdb=10 * log10( power);
powerdb_hw = 10* log10(power_hw);

f5 = figure;
plot( surface_abscurrent(:,1), powerdb(:,1));
title("E Value for Reception"); %Quarter 
ylabel("Electric Field(E)");
xlabel("Distance(X)in Meters");

f6 = figure;
plot( surface_abscurrent_hw(:,1), smoothdata(powerdb_hw(:,1)));
title("E Value for Reception");
ylabel("Electric Field(E)");
xlabel("Distance(X)in Meters");

% hankel function

function H = H02(arg)
    H = besselh(0,2,arg);
end

% 
function E = Ei_Rad(arg)
    E = -H02(arg);
end

function X = X(n, qw)
    X = (n-1) * (qw * 13);
end



function X_hw = X_hw(n, hw)
    X_hw = (n-1) * (hw * 13);
end

function Y = Y(n, y, qw)
    Temp = X(n, qw) / 10.0;
    temp1 = floor(Temp);
    prop = Temp - temp1;
    Y = y(temp1+1) + prop * (y(temp1+2)-y(temp1+1));
end

function Y_hw = Y_hw(n, y, hw)
    Temp_hw = X_hw(n, hw) / 10.0;
    temp1_hw = floor(Temp_hw);
    prop_hw = Temp_hw - temp1_hw;
    Y_hw = y(temp1_hw+1) + prop_hw * (y(temp1_hw+2)-y(temp1_hw+1));
end



% R_source_obs function: the distance between source point and observation
% point.
function RSO = R_source_obs(n, qw, y)
    RSO_SQR = (X(n,qw))^2 + (Y(n, y, qw)+2.4-442)^2;
    
    RSO = sqrt(double(RSO_SQR));
end

function RSO_hw = R_source_obs_hw(n, hw, y_hw)
    RSO_SQR_hw = (X(n,hw))^2 + (Y_hw(n, y_hw, hw)+2.4-442)^2;
    
    RSO_hw = sqrt(double(RSO_SQR_hw));
end

% R_surf_obs function: the distance between surface to observation point
% Params: 1. n -> temp1 number of the surface point
%         2. m -> temp1 number of the observation point
function RSUO = R_surface_obs(n, m, qw, y)
    RSUO_SQR = (X(n, qw) - X(m, qw))^2 + (Y(m, y, qw) +2.4- Y(n, y, qw))^2;
    RSUO = sqrt(double(RSUO_SQR));
end

function RSUO_hw = R_surface_obs_hw(n, m, hw, y_hw)
    RSUO_SQR_hw = (X(n, hw) - X(m, hw))^2 + (Y(m, y_hw, hw) +2.4- Y_hw(n, y_hw, hw))^2;
    RSUO_hw = sqrt(double(RSUO_SQR_hw));
end

% R_p_q function: distance between point
function RPQ = R_p_q(m, n, qw, y)
    RPQ_SQR = (X(m, qw) - X(n, qw))^2  + (Y(m, y, qw) - Y(n, y, qw))^2;
    RPQ = sqrt(RPQ_SQR);
end

% R_p_q function: distance between point
function RPQ_hw = R_p_q_hw(m, n, hw, y_hw)
    RPQ_SQR_hw = (X(m, hw) - X(n, hw))^2  + (Y(m, y_hw, hw) - Y(n, y_hw, hw))^2;
    RPQ_hw = sqrt(RPQ_SQR_hw);
end

% R_soruce_p function: distance between source and the point
function RSP = R_source_p(n, qw, y)
    RSP_SQR = (0 - X(n, qw))^2 + (442 - Y(n, y, qw))^2;
    RSP = sqrt(RSP_SQR);
end

function RSP_hw = R_source_p_hw(n, hw, y_hw)
    RSP_SQR_hw = (0 - X(n, hw))^2 + (442 - Y_hw(n, y_hw, hw))^2;
    RSP_hw = sqrt(RSP_SQR_hw);
end




