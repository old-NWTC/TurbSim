%http://wind.nrel.gov/MetData/135mData/M4Twr/20Hz/mat/2014/01/13/01134_17_40_00_013.mat
% http://wind.nrel.gov/MetData/135mData/M4Twr/20Hz/mat/2014/09/30/09304_05_10_00_030.mat

data_loc  = 'Y:/wind/windweb/MetData/135mData/M4Twr/20Hz/mat';
%data_file = [ data_loc '/2014/01/13/01134_17_40_00_013.mat'];
%input_file = '20140113_174000_inertial.TimeSer';
%heights=[15 30 50 76 100]; %m %131 

% data_file = [ data_loc '/2014/09/29/09294_05_10_00_029.mat'];
% input_file = '20140929_051000_inertial.TimeSer';
% heights=[15 30 50 76 100 131]; %m 


data_file = [ data_loc '/2014/01/13/01134_16_40_00_013.mat'];
input_file = '20140113_164000_inertial_new.TimeSer';
heights=[15 30 50 76 100]; %m 


% data_file = [ data_loc '/2014/01/12/01124_16_40_00_012.mat'];
% input_file = '20140112_164000_inertial.TimeSer';
% heights=[15 30 50 76 100 131]; %m 


load(data_file);

nPoints = length(heights);

% initialize
nt = length(Sonic_dt_clean_15m.val);
t  = zeros(nt,  nPoints);
uvw= zeros(nt,3*nPoints);

RefID = 4;

% in inertial reference frame coordinates:
for iz = 1:length(heights)    
    thisHeight =  num2str(heights(iz));
    t(:,iz) = eval(['Sonic_dt_clean_' thisHeight 'm.val']);
    
    indx = (iz-1)*3+1;
    uvw(:,indx)   = eval(['Sonic_x_clean_' thisHeight 'm.val']); %u
    uvw(:,indx+1) = eval(['Sonic_y_clean_' thisHeight 'm.val']); %v
    uvw(:,indx+2) = eval(['Sonic_z_clean_' thisHeight 'm.val']); %w            
end

y = zeros(size(heights));

%% write input file
writeTurbSimTimeSeriesInput(input_file, data_file, uvw, t, y, heights, RefID);
fclose all;

if any(isnan(uvw(:)))
    error('invalid data in output')
end

%%
% for fun, let's just guess at the coherence function:
[a,b,alpha_h, alpha_v] = InertialData_to_Coherence(uvw,t,y,heights);
%%
fprintf('\n\nEstimated values for primary TurbSim input file:\n');
comp='uvw';
for i=1:3
fprintf( '%s-component coherence parameters for general model are ("%f, %f")\n', comp(i), a(i), b(i) );
end
fprintf( 'CohExp = 0.0\n')
%fprintf( 'HFlowAng = %s degrees\n', -mean(alpha_h) );
%fprintf( 'VFlowAng = %s degrees\n', -mean(alpha_v) );
fprintf( '\n');

%% rotated

    %% let's put this in a different format so we can rotate it easier:
    % probably could just use reshape...
    velocity = zeros(nt,3,1,nPoints);
    IndxStart=0;
    for ip=1:nPoints        
        velocity(:,:,1,ip) = uvw(:,IndxStart+(1:3));
        IndxStart = IndxStart + 3;
    end
    
figure; 
subplot(3,1,1); plot(t,uvw(:,1:3:end)); title('x');
subplot(3,1,2); plot(t,uvw(:,2:3:end)); title('y');
subplot(3,1,3); plot(t,uvw(:,3:3:end)); title('z');

    [velocity, alpha_h2, alpha_v2] = RotateVelocityComponents( velocity, alpha_h(RefID), alpha_v(RefID));
    
    IndxStart=0;
    for ip=1:nPoints        
        uvw(:,IndxStart+(1:3)) = velocity(:,:,1,ip);
        IndxStart = IndxStart + 3;
    end 

figure; 
subplot(3,1,1); plot(t,uvw(:,1:3:end)); title('u');
subplot(3,1,2); plot(t,uvw(:,2:3:end)); title('v');
subplot(3,1,3); plot(t,uvw(:,3:3:end)); title('w');
    
%% write input file
writeTurbSimTimeSeriesInput(['meanDirRemoved' input_file], data_file, uvw, t, y, heights, 4);
fclose all;
    