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
fid = fopen(input_file,'wt');
endline = char([10]); %'\r\n';

fprintf(fid,'--------------TurbSim v2.00.* User Time Series Input File-----------------------%s', endline);
fprintf(fid,'Time series input from %s. Using rotated series %s', data_file, endline);
fprintf(fid,'--------------------------------------------------------------------------------%s', endline);
fprintf(fid,' True        containsW  - Determines whether the time series data includes w-component winds (flag) [if false, w will be simulated using inputs in the main TurbSim input file]%s', endline);
fprintf(fid,' %10.0f  nPoints    - Number of time series points contained in this file (-)%s', length(heights), endline);
fprintf(fid,'PointID   Pointyi     Pointzi     ! PointID is my current thought for a way to specify order of importance%s', endline);
fprintf(fid,'   (-)     (m)         (m)%s', endline);
for iz = 1:length(heights)   
    fprintf( fid, '%5.0f   %10.5f  %10.5f %s', iz, 0.0, heights(iz), endline );   
end
fprintf(fid,'--------Time Series-------------------------------------------------------------%s', endline);
fprintf(fid,' Elapsed Time     ' );
for iz=1:length(heights)
    fprintf (fid, 'Point%02.0f%s_%03.0fm ', iz, 'u', heights(iz), iz, 'v', heights(iz), iz, 'w', heights(iz) );
end
fprintf( fid, '%s', endline);
fprintf(fid,'          (s)     ' );
for iz=1:length(heights)
    for ic=1:3
    fprintf (fid, '    (m/s)     '  );
    end
end
fprintf( fid, '%s', endline);

fmtStr = repmat('%13.4f ',1,size(uvw,2)+1);
fmtStr = [fmtStr '%s'];

for it = 1:size(uvw,1)   
    fprintf( fid, fmtStr, t(it,1), uvw(it,:), endline );   
end
fclose all;