


clc, clear 

tx_max=1;
rx_max=32;

%% Create a batch file 
filename=sprintf('automate1x32.sh');
fid= fopen (filename, 'w');

fprintf(fid, '#!/bin/bash -x \n');
fprintf(fid, '#SBATCH --time=01:00:00 \n');
fprintf(fid, '#SBATCH --mail-user=a.mozaffari@fz-juelich.de \n');
fprintf(fid, '#SBATCH --mail-type=ALL \n');




fprintf(fid, '\n');fprintf(fid, '\n');

for tx_num=1:tx_max
    
    formatSpec= 'cd TX_number_%d \n';
    fprintf(fid,formatSpec,tx_num);
    
    
    for rx_num=1:rx_max
               
    formatSpec= 'cd TX_number_%d_RX_number_%d \n';
    fprintf(fid,formatSpec,tx_num,rx_num);
        
    fprintf(fid,'chmod +x GprMax3D_omp \n');


    formatSpec= 'sbatch sbatch_tx%d_rx%d.sh \n';
    fprintf(fid,formatSpec,tx_num,rx_num);

    fprintf(fid,'cd .. \n');
      
    end    
  fprintf(fid,'cd .. \n');
end
