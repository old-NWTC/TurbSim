template = [ 'C:\Users\bjonkman\Documents\DATA\DesignCodes\preprocessors\TurbSim\' ...
             'SVNdirectory\branches\Modularization\ExampleFiles\TurbSim.inp' ]; 
TestFiles = {'Kaimal.inp'
             'vonKarm.inp'
             'Kaimal_15.inp'
             'vonKarm_15.inp'
             'smooth.inp'
             'KHTest.inp'
             'GPLLJ.inp'
             'HYDRO_TIDAL.inp'
             'GPLLJ_Large.inp'
             'IECKAI_Large.inp'
             'TurbSim_UsrSpec.inp'
             'UsrShear.inp'
             'UsrVkm.inp'             
             'UserTimeSeries.inp' };
             
% TestFiles = {'UserTimeSeries\Kaimal_UserTS.inp'};
% TestFiles = {'..\ExampleFiles\TurbSim_Hydro.inp'};


for i=1:length(TestFiles)

    oldfile = [ 'C:\Users\bjonkman\Documents\DATA\DesignCodes\preprocessors\TurbSim\SVNdirectory\' ...
                'branches\Modularization\CertTest\' TestFiles{i}];
    TSpar = Fast2Matlab(oldfile,2);
    Matlab2FAST(TSpar, template, oldfile, 2); %contains 2 header lines    
end
             