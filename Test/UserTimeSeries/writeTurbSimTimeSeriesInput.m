function writeTurbSimTimeSeriesInput(input_file, data_file, uvw, t, y, z)
% size uvw is (#timeSteps, 3*length(z))

if isempty(y) 
    y = zeros(size(z));
elseif isscalar(y)
    y = ones(size(z))*y;
end
    

fid = fopen(input_file,'wt');
endline = char([10]); %'\r\n';

fprintf(fid,'--------------TurbSim v2.00.* User Time Series Input File-----------------------%s', endline);
fprintf(fid,'Time series input from %s. Using rotated series %s', data_file, endline);
fprintf(fid,'--------------------------------------------------------------------------------%s', endline);
fprintf(fid,' %10.0f  nComp    - Number of velocity components in the file (1=u component only; 2=u & v components; 3=u,v,w) [if < 3 other components will be generated using values from input file]%s', 3, endline);
fprintf(fid,' %10.0f  nPoints  - Number of time series points contained in this file (-)%s', length(z), endline);
fprintf(fid,'PointID   Pointyi     Pointzi     ! PointID is my current thought for a way to specify order of importance%s', endline);
fprintf(fid,'   (-)     (m)         (m)%s', endline);
for iz = 1:length(z)   
    fprintf( fid, '%5.0f   %10.5f  %10.5f %s', iz, y(iz), z(iz), endline );   
end
fprintf(fid,'--------Time Series-------------------------------------------------------------%s', endline);
fprintf(fid,' Elapsed Time     ' );
for iz=1:length(z)
    fprintf (fid, 'Point%02.0f%s_%03.0fm ', iz, 'u', z(iz), iz, 'v', z(iz), iz, 'w', z(iz) );
end
fprintf( fid, '%s', endline);
fprintf(fid,'          (s)     ' );
for iz=1:length(z)
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
fclose(fid);

return
end
