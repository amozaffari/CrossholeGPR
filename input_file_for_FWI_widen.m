clc; clear;


% function that produce the input file for FWI 


%% these data will be read in 

%0.1 _____ FDTD setup 

 domain_size_x=7.02;
 domain_size_y=0.90;
 domain_size_z=11.07;

 mesh_x=0.01; 
 mesh_y=0.01;
 mesh_z=0.01;
 nips_number=38408;
 
 k_factor=9; % ratio of the eps and sig model to the model discritization 
 
 
 n_dx=round(domain_size_x/mesh_x);  n_dy=round(domain_size_y/mesh_y);  n_dz=round(domain_size_z/mesh_z);  % cell numbers in each axis 

 time_window=200e-9;
 pml_number=20;

 %0.2 _____ Crosshole setup info 

% Boreholes configurations
% TX borehole refers as borehole loacted in left and RX in right 
water_table       =  6.84;    % water table  ( where the saturated and unsaturated zone changed)
unsaturated_depth = domain_size_z - water_table;

tx_x_borehole     =  0.979;              % x coordination of the borehole TX -  recip 1 
%tx_x_borehole     =  5.910;              % x coordination of the borehole TX - recip 2 

tx_y_borehole     =  domain_size_y/2;    % y coordination of the borehole TX  
tx_z1_borehole    =  0.2;                % bottom the borehole TX 
tx_z2_borehole    =  domain_size_z - 0.2;     % top of the borehole 
bh_tx_r_1         =  0.06;               % borehole raduis 1 
rx_x_borehole     =  5.910;              % x coordination of the borehole RX - recip 1 
%rx_x_borehole     =  0.979;              % x coordination of the borehole RX - recip 2

rx_y_borehole     =  domain_size_y/2;    % y coordination of the borehole RX 
rx_z1_borehole    =  0.2;                 % top of the borehole RX 
rx_z2_borehole    =  domain_size_z - 0.2;                 % depth of the borehole TX 
bh_rx_r_1         =  0.06;               % borehole raduis 1 

%0.3 _____ Crosshole configuration  

tx_max=11;           %  widen data 22
rx_max=57;           % widen data 114
offset_tx=0.5;      % distance between txs position
offset_rx=0.1;      % distance between rxs position 
first_tx_depth = domain_size_z - 9.656; % depth of the first  measurment _ recip1 
first_rx_depth = domain_size_z - 9.861; % depth of the first  measurment - recip 1

%first_tx_depth = domain_size_z - 9.789; % depth of the first  measurment _ recip2 
%first_rx_depth = domain_size_z - 9.927; % depth of the first  measurment - recip 2  



%0.4 _____ Antennas detil (Sensor and software)
center_frequency=92; % MHz antenna 
ant_length=1.29 ;    % antenna length 
ant_fp=0.26;         % feeding point distance to atenna tip 
ant_half_pec=0.25;   % antenna half PEC length (one-arm)
ant_r1=0.03;         % buddy of the antenna 
ant_r2=0.02;         % resisitors and the wiring
ant_rs_length=0.01;  %length of each antenna segment
ant_rs_offset=0.01;  % distance between two resisitor segments 
ant_fp_length=0.01;  % antenna feeding point vaccum from each side 


%% creating dir 

    folder_name_tx=sprintf('input_FWI');  % for every TX number, it will creat a new folder 
    mkdir (folder_name_tx) 
    cd (folder_name_tx)

%% coordinate file
% -- create the GPRMAX input file 
coordinate_file_name=sprintf('trnrec_FWI.cor');
fid= fopen (coordinate_file_name, 'w');
fprintf(fid, '--------------------- Transmitter-Coordinates ----------------------- \n');
fprintf(fid,'%d',2*tx_max);
fprintf(fid, '\n');  % number of transmitter in recipocal positions 

% recipocal 1 TX position 

first_tx_depth=4.656;

for tx_num=1:tx_max              % TX loop    
tx_position=((tx_num-1)*offset_tx)+first_tx_depth;
formatSpec= '%4.3f  %4.3f \n';  
fprintf(fid,formatSpec,tx_x_borehole, tx_position);
fprintf(fid,'1'); fprintf(fid, '\n');

end 
% Recipocal 2 

first_tx_depth=4.789;

for tx_num=1:tx_max              % TX loop 
tx_position=((tx_num-1)*offset_tx)+first_tx_depth;
formatSpec= '%4.3f  %4.3f \n';  
fprintf(fid,formatSpec,rx_x_borehole, tx_position);
fprintf(fid,'1'); fprintf(fid, '\n');

end 

fprintf(fid, '--------------------- Receiver-Coordinates ----------------------- \n');
fprintf(fid,'%d',2*rx_max);
fprintf(fid, '\n');  % number of transmitter in recipocal positions 

first_rx_depth= 4.278;
for rx_num=1:rx_max          % RX loop 
rx_position=((rx_num-1)*offset_rx)+first_rx_depth;
formatSpec= '%4.3f %4.3f \n';  
fprintf(fid,formatSpec,rx_x_borehole, rx_position);
fprintf(fid,'1'); fprintf(fid, '\n');
        
end 


first_rx_depth=4.344;
for rx_num=1:rx_max          % RX loop 
rx_position=((rx_num-1)*offset_rx)+first_rx_depth;
formatSpec= ' %4.3f %4.3f \n';  
fprintf(fid,formatSpec,tx_x_borehole, rx_position);
fprintf(fid,'1'); fprintf(fid, '\n');
        
end 



%% Mark.E file 

makE_file_name=sprintf('mark_FWI.E');
fid1= fopen (makE_file_name, 'w');
fprintf(fid1, 'Include(=1.0)/Exclude(=0.0) Traces \n');

for tx_num=1:tx_max              % Recipocal 1     
    formatSpec= 'TNr.:%1d \n';  
    fprintf(fid1,formatSpec,tx_num);   
    for rx_num=1:rx_max          % RX loop with included traces 
    fprintf(fid1,'1'); fprintf(fid1, '\n');        
    end    
    for rx_num=1:rx_max          % RX loop with included traces 
    fprintf(fid1,'0'); fprintf(fid1, '\n');        
    end             
end 

for tx_num=tx_max+1:2*tx_max              % Recipocal 1 
    
    formatSpec= 'TNr.:%1d \n';  
    fprintf(fid1,formatSpec,tx_num);    
    for rx_num=1:rx_max          % RX loop with included traces 
    fprintf(fid1,'0'); fprintf(fid1, '\n');        
    end     
    for rx_num=1:rx_max          % RX loop with included traces 
    fprintf(fid1,'1'); fprintf(fid1, '\n');        
    end            
end 





cd ..