 clc; clear;close all;

 %%  input 
convert_ms_mns=1e-9;   % convert the ms-1 to the mns-1
 
 
 
 %%
 
% Change 

%Desktop
%direct='/Users/a.mozaffari/Dropbox/IBG-3/Results/Antenna_model/0_FWI_GPRMAX/synthtic_model';
% MacBook
%direct='//Users/amirpashamozaffari/Dropbox/IBG-3/Results/Antenna_model/0_FWI_GPRMAX/synthtic_model';

tx_max=1;
rx_max=57;
TX_plot=1;

samples = 10386; % read from the header file ---> iteration

Field.ex=zeros(samples,tx_max,rx_max);
Field.ey=zeros(samples,tx_max,rx_max);
Field.ez=zeros(samples,tx_max,rx_max);
Field.t=zeros(samples,tx_max,rx_max);

% NOTE : Directory need to be changed manually in gprMax copy part 
%cd (direct) 
 


for tx_num=1:tx_max              % TX loop 

    
    % -- go to direcxtoy TX

%     folder_name_tx=sprintf('TX_number_%d',tx_num);      
%     cd (folder_name_tx)
    
    variable_name=sprintf('Fields%d.x',tx_num);
    
    for rx_num=1:rx_max          % RX loop 
    
                % -- go to directory RX  
                folder_name_rx=sprintf('TX_number_%d_RX_number_%d',tx_num,rx_num);  
                cd (folder_name_rx)

                
%% Extract the information of the recorded trace 
      
                % -- open the GPRMAX input file 
                gprmax_input_name=sprintf('tx%d_rx%d.out',tx_num,rx_num);
                % fid= fopen (gprmax_input_name, 'r')
                [Header1,Field1]=gprmax(gprmax_input_name);
                
                
                for sample_num=1:samples
                                
                 Field.ex(sample_num,tx_num,rx_num)=Field1.ex(sample_num,1);
                 Field.ey(sample_num,tx_num,rx_num)=Field1.ey(sample_num,1);
                 Field.ez(sample_num,tx_num,rx_num)=Field1.ez(sample_num,1);
                 Field.t(sample_num,tx_num,rx_num)=Field1.t(sample_num,1);  
                 
                end 
                
                cd ..
    end 

cd ..

end 
%% Save 
%save('All_data','Field');

%% Creat the shot-gather 

% Pre-alocate 
fprintf('Alloctaing the MTX \n');

Shot_gather_ex=zeros(samples,rx_max);
Shot_gather_ey=zeros(samples,rx_max);
Shot_gather_ez=zeros(samples,rx_max);
Shot_gather_t=zeros(samples,rx_max);


Shot_gather_norm_ex=zeros(samples,rx_max);
Shot_gather_norm_ey=zeros(samples,rx_max);
Shot_gather_norm_ez=zeros(samples,rx_max);

Image_scan_ex=zeros(samples,rx_max);
Image_scan_ey=zeros(samples,rx_max);
Image_scan_ez=zeros(samples,rx_max);

Image_scan_norm_ex=zeros(samples,rx_max);
Image_scan_norm_ey=zeros(samples,rx_max);
Image_scan_norm_ez=zeros(samples,rx_max);

% Shot - gather 


k_mag=1.0;

for  rx_num=1:rx_max 

        for sample_num=1:samples
            
        % absolute shot gather     
        Shot_gather_ex(sample_num,rx_num)=Field.ex(sample_num,TX_plot,rx_num)+rx_num;
        Shot_gather_ey(sample_num,rx_num)=Field.ey(sample_num,TX_plot,rx_num)+rx_num;
        Shot_gather_ez(sample_num,rx_num)=Field.ez(sample_num,TX_plot,rx_num)+rx_num;
        Shot_gather_t(sample_num,rx_num)=Field.t(sample_num,TX_plot,rx_num)+rx_num;
       
        % normalized shot gather 
        
        Shot_gather_norm_ex(sample_num,rx_num)=(Field.ex(sample_num,TX_plot,rx_num)/max(abs(Field.ex(:,TX_plot,rx_num)))*k_mag)+rx_num;
        Shot_gather_norm_ey(sample_num,rx_num)=(Field.ey(sample_num,TX_plot,rx_num)/max(abs(Field.ey(:,TX_plot,rx_num)))*k_mag)+rx_num;
        Shot_gather_norm_ez(sample_num,rx_num)=(Field.ez(sample_num,TX_plot,rx_num)/max(abs(Field.ez(:,TX_plot,rx_num)))*k_mag)+rx_num;        
        
        
        
        %% todo list tomorrow : 
        
        
        % Absolute image scan 
        
        Image_scan_ex(sample_num,rx_num)=Field.ex(sample_num,TX_plot,rx_num);
        Image_scan_ey(sample_num,rx_num)=Field.ey(sample_num,TX_plot,rx_num);
        Image_scan_ez(sample_num,rx_num)=Field.ez(sample_num,TX_plot,rx_num);
        
        % normalized image scan 
        
        Image_scan_norm_ex(sample_num,rx_num)=(Field.ex(sample_num,TX_plot,rx_num)/max(abs(Field.ex(:,TX_plot,rx_num)))*k_mag);
        Image_scan_norm_ey(sample_num,rx_num)=(Field.ey(sample_num,TX_plot,rx_num)/max(abs(Field.ey(:,TX_plot,rx_num)))*k_mag);
        Image_scan_norm_ez(sample_num,rx_num)=(Field.ez(sample_num,TX_plot,rx_num)/max(abs(Field.ez(:,TX_plot,rx_num)))*k_mag);        
        
        
        
        % frequncy spectrum of the traces 
        
        
        
        
        
        
        end
        
        formatSpec= 'Extraction and normalization of the TX: %d and RX: %d trace is done\n'; % antenna buddy 
        fprintf(formatSpec,TX_plot,rx_num);

end




%% Plot

%% Absolute 

% EX 
ex=figure; 
    plot(Field1.t,Shot_gather_ex,...
        'Color','black',...
        'LineWidth',1);
    
    
    grid on

    xhandle=xlabel('X axis');
    yhandle=ylabel('Y axis');
    set(ex, 'Position', [1 1 800 900])
    set(xhandle,'FontSize',17)
    set(yhandle,'FontSize',17)
    xlim([0 2E-07]);
    ylim([0 46])
    ylabel ('RX number');
    zlabel ('Time');    
    title(['Absolute amplitude shot gather of transmitter number:',num2str(TX_plot),'for Ex field'],'FontSize',14) 
    fig_title=sprintf('EX_shot_gather_TX_number_%d',TX_plot);  

    savefig(ex,fig_title); 
  %  saveas(ex,fig_title,'jpeg');    
    
  
 % EY 
ey=figure; 
    plot(Field1.t,Shot_gather_ey,...
        'Color','black',...
        'LineWidth',1);
    
    
    grid on

    xhandle=xlabel('X axis');
    yhandle=ylabel('Y axis');
    set(ey, 'Position', [1 1 800 900])
    set(xhandle,'FontSize',17)
    set(yhandle,'FontSize',17)
    xlim([0 2E-07]);
    ylim([0 46])
    ylabel ('RX number');
    zlabel ('Time');    
    title(['Absolute amplitude shot gather of transmitter number:',num2str(TX_plot),'for Ey field'],'FontSize',14) 
    fig_title=sprintf('EY_shot_gather_TX_number_%d',TX_plot);  

    savefig(ey,fig_title); 
    % saveas(ey,fig_title,'jpeg');   
    
    % Ez 
ez=figure; 
    plot(Field1.t,Shot_gather_ez,...
        'Color','black',...
        'LineWidth',1);
    
    
    grid on

    xhandle=xlabel('X axis');
    yhandle=ylabel('Y axis');
    set(ez, 'Position', [1 1 800 900])
    set(xhandle,'FontSize',17)
    set(yhandle,'FontSize',17)
    xlim([0 2E-07]);
    ylim([0 46])
    ylabel ('RX number');
    zlabel ('Time');    
    title(['Absolute amplitude shot gather of transmitter number:',num2str(TX_plot),'for Ez field'],'FontSize',14) 
    fig_title=sprintf('Ez_shot_gather_TX_number_%d',TX_plot);  

    savefig(ez,fig_title); 
   % saveas(ez,fig_title,'jpeg');  


%% Normalized 
% EX 
ex=figure; 
    plot(Field1.t,Shot_gather_norm_ex,...
        'Color','black',...
        'LineWidth',1);
      
    grid on

    xhandle=xlabel('X axis');
    yhandle=ylabel('Y axis');
    set(ex, 'Position', [1 1 800 900])
    set(xhandle,'FontSize',17)
    set(yhandle,'FontSize',17)
    xlim([0 2E-07]);
    ylim([0 46])
    ylabel ('RX number');
    zlabel ('Time');    
    title(['Normalized amplitude shot gather of transmitter number:',num2str(TX_plot),'for Ex field'],'FontSize',14) 
    fig_title=sprintf('EX_normalized_shot_gather_TX_number_%d',TX_plot);  

    savefig(ex,fig_title); 
 %   saveas(ex,fig_title,'jpeg');    
    
  
 % EY 
ey=figure; 
    plot(Field1.t,Shot_gather_norm_ey,...
        'Color','black',...
        'LineWidth',1);
    
    
    grid on

    xhandle=xlabel('X axis');
    yhandle=ylabel('Y axis');
    set(ey, 'Position', [1 1 800 900])
    set(xhandle,'FontSize',17)
    set(yhandle,'FontSize',17)
    xlim([0 2E-07]);
    ylim([0 46])
    ylabel ('RX number');
    zlabel ('Time');    
    title(['Normalized amplitude shot gather of transmitter number:',num2str(TX_plot),'for Ey field'],'FontSize',14) 
    fig_title=sprintf('EY_normalized_shot_gather_TX_number_%d',TX_plot);  

    savefig(ey,fig_title); 
  %  saveas(ey,fig_title,'jpeg');   
    
    % Ez 
ez=figure; 
    plot(Field1.t,Shot_gather_norm_ez,...
        'Color','black',...
        'LineWidth',1);
    
    
    grid on

    xhandle=xlabel('X axis');
    yhandle=ylabel('Y axis');
    set(ez, 'Position', [1 1 800 900])
    set(xhandle,'FontSize',17)
    set(yhandle,'FontSize',17)
    xlim([0 2E-07]);
    ylim([0 46])
    ylabel ('RX number');
    zlabel ('Time');    
    title(['Normalized amplitude shot gather of transmitter number:',num2str(TX_plot),'for Ez field'],'FontSize',14) 
    fig_title=sprintf('Ez_normalized_shot_gather_TX_number_%d',TX_plot);  

    savefig(ez,fig_title); 
 %   saveas(ez,fig_title,'jpeg');   
    
    
  %% image scan   

  % Ex 
ex=figure; 
    imagesc(flipud(Image_scan_ex'));
    title(['Absolute amplitude gather of transmitter number:',num2str(TX_plot),'for Ex field'],'FontSize',14) 
    fig_title=sprintf('Ex_scan_gather_TX_number_%d',TX_plot);        
    colorbar;
    ylabel('RX number') ;xlabel('Sample') ;    
    hold on
    colormap ('Jet')
    savefig(ex,fig_title); 
%    saveas(ex,fig_title,'jpeg');
    
 % Ey 
ey=figure; 
    imagesc(flipud(Image_scan_ey'));
    title(['Absolute amplitude gather of transmitter number:',num2str(TX_plot),'for Ey field'],'FontSize',14) 
    fig_title=sprintf('Ey_scan_gather_TX_number_%d',TX_plot);        
    colorbar;
    ylabel('RX number') ;xlabel('Sample') ;    
    hold on
    colormap ('Jet')
    savefig(ey,fig_title); 
 %   saveas(ey,fig_title,'jpeg'); 
    
% Ez 
ez=figure; 
    imagesc(flipud(Image_scan_ez'));
    title(['Absolute amplitude gather of transmitter number:',num2str(TX_plot),'for Ez field'],'FontSize',14) 
    fig_title=sprintf('Ez_scan_gather_TX_number_%d',TX_plot);        
    colorbar;
    ylabel('RX number') ;xlabel('Sample') ;    
    hold on
    colormap ('Jet')   
    savefig(ez,fig_title); 
 %   saveas(ez,fig_title,'jpeg');     
    
% Normalized 

  % Ex 
ex=figure; 
    imagesc(flipud(Image_scan_norm_ex'));
    title(['Normalized amplitude gather of transmitter number:',num2str(TX_plot),'for Ex field'],'FontSize',14) 
    fig_title=sprintf('Ex_norm_scan_gather_TX_number_%d',TX_plot);        
    colorbar;
    ylabel('RX number') ;xlabel('Sample') ;    
    hold on
    colormap ('Jet')
    savefig(ex,fig_title); 
    saveas(ex,fig_title,'jpeg');
    
 % Ey 
ey=figure; 
    imagesc(flipud(Image_scan_norm_ey'));
    title(['Normalized amplitude gather of transmitter number:',num2str(TX_plot),'for Ey field'],'FontSize',14) 
    fig_title=sprintf('Ey_norm_scan_gather_TX_number_%d',TX_plot);        
    colorbar;
    ylabel('RX number') ;xlabel('Sample') ;    
    hold on
    colormap ('Jet')
    savefig(ey,fig_title); 
 %   saveas(ey,fig_title,'jpeg'); 
    
% Ez 
ez=figure; 
    imagesc(flipud(Image_scan_norm_ez'));
    title(['Normalized amplitude gather of transmitter number:',num2str(TX_plot),'for Ez field'],'FontSize',14) 
    fig_title=sprintf('Ez_norm_scan_gather_TX_number_%d',TX_plot);        
    colorbar;
    ylabel('RX number') ;xlabel('Sample') ;    
    hold on
    colormap ('Jet')   
    savefig(ez,fig_title); 
%    saveas(ez,fig_title,'jpeg');    

    

%% 

