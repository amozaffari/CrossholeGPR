%function generateTraces(tx,rx) % 
clc, clear 

%% ---------------------- Directory Change ----------------------- %%  

%Desktop
%direct='/Users/a.mozaffari/Dropbox/IBG-3/Results/Antenna_model/0_FWI_GPRMAX/synthtic_model';
% MacBook
direct='//Users/amirpashamozaffari/Dropbox/IBG-3/Results/Antenna_model/0_FWI_GPRMAX/synthtic_model/3D_cube_build';
cd (direct)

%% ---------------------- TO KNOW ------------------------------ %% 

 E0 = 8.85418781762039080*1e-12;  % convert the abs eps to relative eps 

%% outline of the process 

% read the starting model eps and sigma from FWI
% assign material of each cell a code 
% build the each cube material in the textfile 


% combine it with TX and RX builder from FWI_input generator code 

%% -------------------------- 0.intial info -------------------------- %%

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

%tx_x_borehole     =  0.979;              % x coordination of the borehole TX -  recip 1 
tx_x_borehole     =  5.910;              % x coordination of the borehole TX - recip 2 

tx_y_borehole     =  domain_size_y/2;    % y coordination of the borehole TX  
tx_z1_borehole    =  0.2;                % bottom the borehole TX 
tx_z2_borehole    =  domain_size_z - 0.2;     % top of the borehole 
bh_tx_r_1         =  0.06;               % borehole raduis 1 
%rx_x_borehole     =  5.910;              % x coordination of the borehole RX - recip 1 
rx_x_borehole     =  0.979;              % x coordination of the borehole RX - recip 2

rx_y_borehole     =  domain_size_y/2;    % y coordination of the borehole RX 
rx_z1_borehole    =  0.2;                 % top of the borehole RX 
rx_z2_borehole    =  domain_size_z - 0.2;                 % depth of the borehole TX 
bh_rx_r_1         =  0.06;               % borehole raduis 1 

%0.3 _____ Crosshole configuration  

tx_max=11;           %  widen data 22
rx_max=57;           % widen data 114
offset_tx=0.5;      % distance between txs position
offset_rx=0.1;      % distance between rxs position 


%first_tx_depth = domain_size_z - 9.656; % depth of the first  measurment _ recip1 
%first_rx_depth = domain_size_z - 9.861; % depth of the first  measurment - recip 1

first_tx_depth = domain_size_z - 9.789; % depth of the first  measurment _ recip2 
first_rx_depth = domain_size_z - 9.927; % depth of the first  measurment - recip 2  



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


%% ----------------------- 1. Building up the 3D Cube -------------------- %%

% 1.1 Read the EPS and Sigma 
% caution unit transformation

% 1.1.EPS 
Fid1= fopen('model_synth.eps','rb');
nx_1=fread(Fid1,1,'int');
nz_1=fread(Fid1,1,'int');
EV=fread(Fid1,nx_1*nz_1,'double');fclose(Fid1);
EV=reshape(EV,nz_1,nx_1);
EV=EV/E0;  % realative EPS. 
EV_original=EV;
EV=flipud(EV); % flipped to build the cube from botoom 

% 1.1.2 SIG

Fid2= fopen('model_synth.sig','rb');
nx_1=fread(Fid2,1,'int');
nz_1=fread(Fid2,1,'int');
SA=fread(Fid2,nx_1*nz_1,'double');fclose(Fid2);
SA=reshape(SA,nz_1,nx_1);
SA_original=SA;
SA=flipud(SA);  % flipped to build the cube from botoom 
 
fclose('all');

% 1.2 Buildup the media MTX                    
% 1.2.1 initiation of media_file  

media.f1=zeros(nz_1,nx_1);  % relative permitivity (Epsilon_r)
media.f2=zeros(nz_1,nx_1);  % relative permitivity in infinite frequency (Epsilon_r_inif)
media.f3=zeros(nz_1,nx_1);  % relaxition time in seconds 
media.f4=zeros(nz_1,nx_1);  % conductivity in Siemens/Meters (Sigma)
media.f5=zeros(nz_1,nx_1);  % relative permieability (miue_r)
media.f6=zeros(nz_1,nx_1);  % magnetic cnductivity (Sigma_star) 
media.str=zeros(nz_1,nx_1); % meduim_id 


%1.2.2 filling up the media file 
med_id_max= nz_1*nx_1 + 8;%size(media.f1,1) *size(media.f1,2)  + 10; % #test
med_id=1;  % #test

for z=1:size(media.f1,1) 
    for x=1:size(media.f1,2)            
      media.f1(z,x) = EV(z,x);
      media.f2(z,x) = 0;
      media.f3(z,x) = 0;      
      media.f4(z,x) = SA(z,x);
      media.f5(z,x) = 1.00;      
      media.f6(z,x) = 0;
      media.str(z,x) = med_id;
      med_id=med_id+1;
    end 
    
end

%% ______________________ 2. write down the input file (test) _________ %%
% 2.1. Creating folder and intial data 

for tx_num=10:tx_max              % TX loop 

  
% 2.1.1 creat TX the directory 

    folder_name_tx=sprintf('TX_number_%d',tx_num);  % for every TX number, it will creat a new folder 
    mkdir (folder_name_tx) 
    
    cd (folder_name_tx)

% 2.1.2 Creating RX directory     
    for rx_num=1:rx_max          % RX loop 

                % RX
                folder_name_rx=sprintf('TX_number_%d_RX_number_%d',tx_num,rx_num);  % for every TX number, it will create a new folder 
                mkdir (folder_name_rx)  
                cd (folder_name_rx)

                
                % gprMAX input file 
                gprmax_input_name=sprintf('input_tx%d_rx%d.in',tx_num,rx_num);
                fid= fopen (gprmax_input_name, 'w');
                
                
                % sbatch file 
                sbatch_input_name=sprintf('sbatch_tx%d_rx%d.sh',tx_num,rx_num);
                fid1= fopen (sbatch_input_name, 'w');                    
    
               
% 2.1.3 calculation of the feeding point position 
tx_position=((tx_num-1)*offset_tx)+first_tx_depth;
rx_position=((rx_num-1)*offset_rx)+first_rx_depth;
                
                           
% 2.1. media decleration


% -- create the GPRMAX input file
gprmax_input_name=sprintf('input_tx%d_rx%d.in',tx_num,rx_num);
fid= fopen (gprmax_input_name, 'w');


% header of the 
fprintf(fid, '\n');fprintf(fid, '\n');
fprintf(fid, '%%____ input file for gprMAX3D V.2  results generated with input_gen_V002_____%% \n');
fprintf(fid, '\n');fprintf(fid, '\n');

% Meduim properties decleration 
fprintf(fid, '------------ Meduim Decleration ------------ \n');
fprintf(fid, '           DC (static) relative permittivity of\n');
fprintf(fid, '           |	    relative permittivity at theoretically infinite frequency\n');
fprintf(fid, '           |	    |		relaxation time of the medium\n');
fprintf(fid, '           |	    |		|		DC (static) conductivity of the medium ?? (Siemens/metre)\n');
fprintf(fid, '           |	    | 		|		|		relative permeability of the medium\n');
fprintf(fid, '           |	    |		|		|		|	    magnetic conductivity of the medium\n');
fprintf(fid, '           |	    |		|		|		|	    |		name of the material \n');

fprintf(fid, '#medium:  1.00	0.000 	0.000 	0.0010 		1.000  0.000 air \n');
fprintf(fid, '#medium:  80.00	0.000 	0.000 	0.05 		1.000  0.000 fresh_water \n');
fprintf(fid, '#medium:  04.0	0.000	0.000	0.00001		1.000  0.000 C_1 \n');
fprintf(fid, '#medium:  04.0	0.000	0.000	0.0 	    1.000  0.000 S_1 \n');
fprintf(fid, '#medium:  3.90    0.000	0.000	1e-14 	    1.000  0.000 PVC \n');


for z=1:size(media.f1,1)  
    for x=1:size(media.f1,2) 
        
        formatSpec= '#medium:   %4.4f   %4.2f   %4.2f   %4.4f   %4.2f   %4.2f  media_%d\n'; % media_MTX 
        fprintf(fid,formatSpec,media.f1(z,x),media.f2(z,x),media.f3(z,x),media.f4(z,x),media.f5(z,x),media.f6(z,x),media.str(z,x));
        
    end     
end


fprintf(fid, '\n');fprintf(fid, '\n');
fprintf(fid, '------------ Setup confuguration ------------ \n');

formatSpec= '#domain:    %4.2f   %4.2f   %4.2f \n'; % x domain size  
fprintf(fid,formatSpec,domain_size_x,domain_size_y,domain_size_z);

formatSpec= '#dx_dy_dz:  %4.2f   %4.2f   %4.2f \n'; % discritization   
fprintf(fid,formatSpec,mesh_x,mesh_y,mesh_z);

formatSpec= '#time_window:  %4.9f\n'; % time window   
fprintf(fid,formatSpec,time_window);

fprintf(fid, '#abc_type: pml \n');
formatSpec= '#pml_layers:  %4.2f \n'; % PML layers   
fprintf(fid,formatSpec,pml_number);

fprintf(fid, '#num_of_procs: 9 \n');

formatSpec= '#nips_number: %d\n';     % nips number   
fprintf(fid,formatSpec,nips_number);

formatSpec= '#number_of_media: %d\n'; % media number  
fprintf(fid,formatSpec,med_id_max);
     
fprintf(fid, '\n');fprintf(fid, '\n');

%% ______________________ 3. subsurface built (3D) _________ %%

fprintf(fid, '------------ Subsurface model cube Widen ------------ \n');

% we build the the blocks from the bottom to the top of the cube 
%______________________ MEDUIM model transfer  _________________________ %

x_position=0.0;
y_position=0.0;
z_position=0.0;

for z=1:size(media.f1,1)  

    for x=1:size(media.f1,2)          
                       
        formatSpec= '#box:   %4.4f   %4.4f   %4.4f   %4.4f   %4.2f   %4.2f  media_%d\n'; % media_MTX 
        fprintf(fid,formatSpec,x_position,y_position,z_position,x_position + k_factor*mesh_x,domain_size_y,z_position + k_factor*mesh_z,media.str(z,x));
        x_position=x_position+k_factor*mesh_x;               
    end     
x_position=0.0;     % reseting the x_position of the cubic 
z_position=z_position+k_factor*mesh_z;     
end

fprintf(fid, '\n');fprintf(fid, '\n');
%______________________ BOREHOLE BUILD _________________________ %
fprintf(fid, '------------ Boreholes 1 (TX)  ------------ \n');fprintf(fid, '\n');
% borehole 1 (TX)
% Water-filling 
formatSpec= '#cylinder: z   %4.3f   %4.3f   %4.3f   %4.3f   %4.3f    fresh_water \n'; % Water-filled (Saturated) part of borehole
fprintf(fid,formatSpec,tx_z1_borehole,water_table,tx_x_borehole,tx_y_borehole,bh_tx_r_1);
% air-filling 
formatSpec= '#cylinder: z   %4.3f   %4.3f   %4.3f   %4.3f   %4.3f    air \n'; % air-filled (Saturated) part of borehole
fprintf(fid,formatSpec,water_table,tx_z2_borehole,tx_x_borehole,tx_y_borehole,bh_tx_r_1);
fprintf(fid, '\n');


fprintf(fid, '------------ Boreholes 2 (RX)  ------------ \n');fprintf(fid, '\n');
% borehole 2 (RX)
% Water-filling 
formatSpec= '#cylinder: z   %4.3f   %4.3f   %4.3f   %4.3f   %4.3f    fresh_water \n'; % Water-filled (Saturated) part of borehole
fprintf(fid,formatSpec,rx_z1_borehole,water_table,rx_x_borehole,rx_y_borehole,bh_rx_r_1);
% air-filling 
formatSpec= '#cylinder: z   %4.3f   %4.3f   %4.3f   %4.3f   %4.3f    air \n'; % air-filled (Saturated) part of borehole
fprintf(fid,formatSpec,water_table,rx_z2_borehole,rx_x_borehole,rx_y_borehole,bh_rx_r_1);
fprintf(fid, '\n');

%% Antenna setup 

fprintf(fid, '------------ RLA Antennas  TX & RX ------------ \n');
fprintf(fid, '\n');

%______________________________ TX Antenna build ___________________ 
fprintf(fid, '--- TX Antenna --- \n');

fprintf(fid, '#voltage_source: 0.1 0.92e8 ricker 0 MyDipole \n'); %NOTE: here is the source that we use in general
fprintf(fid, '\n');

fprintf(fid, 'Antenna isolation: \n');
formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	S_1 \n'; % antenna buddy 
fprintf(fid,formatSpec,(tx_position-ant_fp),(tx_position-ant_fp+ant_length),tx_x_borehole,tx_y_borehole,ant_r1);
fprintf(fid, '\n');
fprintf(fid, 'Antenna wire: \n');
formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	pec \n'; % antenna wire 
fprintf(fid,formatSpec,(tx_position-ant_half_pec),(tx_position+ant_half_pec),tx_x_borehole,tx_y_borehole,ant_r2);
fprintf(fid, '\n');
fprintf(fid, 'Resistors: \n');
for res_count=1:10     
    formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	C_1 \n'; % resistors segments in lower part 
    fprintf(fid,formatSpec,((tx_position-ant_half_pec)+(res_count-1)*(ant_rs_length+ant_rs_offset)),((tx_position-ant_half_pec)+(ant_rs_length)+(res_count-1)*(ant_rs_length+ant_rs_offset)),tx_x_borehole,tx_y_borehole,ant_r2);   
end 
fprintf(fid, '\n');
for res_count=1:10     
    formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	C_1 \n'; % resistors segments in lower part 
    fprintf(fid,formatSpec,((tx_position+ant_half_pec)-(ant_rs_length)-(res_count-1)*(ant_rs_length+ant_rs_offset)),((tx_position+ant_half_pec)-(res_count-1)*(ant_rs_length+ant_rs_offset)),tx_x_borehole,tx_y_borehole,ant_r2);   
end 
fprintf(fid, '\n');
fprintf(fid, 'Antenna feeding point void: \n');
formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	free_space \n'; % antenna void  
fprintf(fid,formatSpec,(tx_position-ant_fp_length),(tx_position+ant_fp_length),tx_x_borehole,tx_y_borehole,ant_r2);
fprintf(fid, '\n');

%___________________________ RX Antenna build ___________________ 

fprintf(fid, '\n');
fprintf(fid, '--- RX Antenna --- \n');
fprintf(fid, 'Antenna isolation: \n');
formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	S_1 \n'; % antenna buddy 
fprintf(fid,formatSpec,(rx_position-ant_fp),(rx_position-ant_fp+ant_length),rx_x_borehole,rx_y_borehole,ant_r1);
fprintf(fid, '\n');
fprintf(fid, 'Antenna wire: \n');
formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	pec \n'; % antenna wire 
fprintf(fid,formatSpec,(rx_position-ant_half_pec),(rx_position+ant_half_pec),rx_x_borehole,rx_y_borehole,ant_r2);
fprintf(fid, '\n');
fprintf(fid, 'Resistors: \n');
for res_count=1:10     
    formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	C_1 \n'; % resistors segments in lower part 
    fprintf(fid,formatSpec,((rx_position-ant_half_pec)+(res_count-1)*(ant_rs_length+ant_rs_offset)),((rx_position-ant_half_pec)+(ant_rs_length)+(res_count-1)*(ant_rs_length+ant_rs_offset)),rx_x_borehole,rx_y_borehole,ant_r2);   
end 
fprintf(fid, '\n');
for res_count=1:10     
    formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	C_1 \n'; % resistors segments in lower part 
    fprintf(fid,formatSpec,((rx_position+ant_half_pec)-(ant_rs_length)-(res_count-1)*(ant_rs_length+ant_rs_offset)),((rx_position+ant_half_pec)-(res_count-1)*(ant_rs_length+ant_rs_offset)),rx_x_borehole,rx_y_borehole,ant_r2);    
end 
fprintf(fid, '\n');
fprintf(fid, 'Antenna feeding point void: \n');
formatSpec= '#cylinder: z  %4.3f   %4.3f   %4.3f  %4.3f    %4.3f	free_space \n'; % feeding point RX 
fprintf(fid,formatSpec,(rx_position-ant_fp_length),(rx_position+ant_fp_length),rx_x_borehole,rx_y_borehole,ant_r2);
fprintf(fid, '\n'); fprintf(fid, '\n');

% starting the simulation  in GPRMAX

formatSpec= '#analysis: 1  tx%d_rx%d.out b \n';
fprintf(fid,formatSpec,tx_num,rx_num );

% TX position
fprintf(fid, '\n');
fprintf(fid, 'Source position: \n');
formatSpec= '#tx: z   %4.3f   %4.3f   %4.3f  MyDipole 0.0 	200e-9 \n';  
fprintf(fid,formatSpec,tx_x_borehole,tx_y_borehole, tx_position);
fprintf(fid, '\n');

% RX position 
fprintf(fid, '\n');
fprintf(fid, 'Sink position: \n');
formatSpec= '#rx:  %4.3f   %4.3f   %4.3f \n';  
fprintf(fid,formatSpec,rx_x_borehole,rx_y_borehole, rx_position);
fprintf(fid, '\n');

fprintf(fid, '#end_analysis:\n');
fprintf(fid, '#messages: y \n');

%% write the Sbatch 

fprintf(fid1, '#!/bin/bash \n');
fprintf(fid1, '#SBATCH --nodes=2 \n');
fprintf(fid1, '#SBATCH --ntasks=2 \n');
fprintf(fid1, '#SBATCH --ntasks-per-node=1 \n');
fprintf(fid1, '#SBATCH --cpus-per-task=9 \n');
fprintf(fid1, '#SBATCH --output=mpi-out.%s \n','%j');
fprintf(fid1, '#SBATCH --error=mpi-err.%s \n','%j');
fprintf(fid1, '#SBATCH --time=06:00:00 \n');
fprintf(fid1, '#SBATCH --partition=batch \n');
fprintf(fid1, '#SBATCH --mail-type=ALL \n');
fprintf(fid1, 'module purge \n');
fprintf(fid1, 'module load GCC \n');
fprintf(fid1, 'module load ParaStationMPI \n');
fprintf(fid1, 'export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \n');

formatSpec= 'srun ./GprMax3D_omp ./input_tx%d_rx%d.in \n';
fprintf(fid1,formatSpec,tx_num,rx_num );

fclose(fid1);

%% change the directory here 
fclose(fid);
    cd ..
    end 
  cd ..  % chnage directory to upper level 
 
  formatSpec= 'Input files are generated for TX: %d and RX: %d \n'; % antenna buddy 
  fprintf(formatSpec,tx_num,rx_num);
    
end 

%% copy the gprmax executable 

for tx_num=1:tx_max              % TX loop   
    for rx_num=1:rx_max          % RX Loop

  %root_name =sprintf('/Users/a.mozaffari/Dropbox/IBG-3/Results/Antenna_model/0_FWI_GPRMAX/synthtic_model/TX_number_%d/TX_number_%d_RX_number_%d',tx_num,tx_num,rx_num);
  root_name =sprintf('/Users/amirpashamozaffari/Dropbox/IBG-3/Results/Antenna_model/0_FWI_GPRMAX/synthtic_model/3D_cube_build/TX_number_%d/TX_number_%d_RX_number_%d',tx_num,tx_num,rx_num);

  copyfile ('GprMax3D_omp',root_name) 

    end 
end




