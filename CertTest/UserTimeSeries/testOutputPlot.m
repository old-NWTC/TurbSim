time_elapsed = (time_UTC.val-time_UTC.val(1) ) *60*60*24;
figure
plot (time_elapsed , Sonic_x_100.val , 'ko' )
hold on
plot (Sonic_dt_clean_100m.val , Sonic_x_clean_100m.val , 'b+' )
plot (Sonic_dt_rotated_100m.val , Sonic_u_100m.val , 'rx' )
legend ( 'Raw data' , ' Cleaned data' , 'Rotated data' )
set ( legend , 'location' , 'best' )