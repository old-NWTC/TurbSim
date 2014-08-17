%http://wind.nrel.gov/MetData/135mData/M4Twr/20Hz/mat/2014/01/13/01134_17_40_00_013.mat

data_loc  = 'Y:/wind/windweb/MetData/135mData/M4Twr/20Hz/mat';
data_file = [ data_loc '/2014/01/13/01134_17_40_00_013.mat'];
input_file = '20140113_174000.inp';



load(data_file);

heights=[76 100 50 30 15]; %m %131 


% initialize
t  = zeros(length(Sonic_dt_rotated_15m.val),length(heights));
uvw= zeros(length(Sonic_dt_rotated_15m.val),3*length(heights));

% in rotated coordinates (along mean wind):
for iz = 1:length(heights)    
    thisHeight =  num2str(heights(iz));
    t(:,iz) = eval(['Sonic_dt_rotated_' thisHeight 'm.val']);
    
    indx = (iz-1)*3+1;
    uvw(:,indx)   = eval(['Sonic_u_' thisHeight 'm.val']); %u
    uvw(:,indx+1) = eval(['Sonic_v_' thisHeight 'm.val']); %v
    uvw(:,indx+2) = eval(['Sonic_w_' thisHeight 'm.val']); %w            
end

%write input file
writeTurbSimTimeSeriesInput(input_file, data_file, uvw, t, y, z);
fclose all;