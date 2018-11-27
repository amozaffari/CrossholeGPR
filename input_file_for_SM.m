
clc,clear;
%% Setup-detail 

 domain_size_x=6;
 domain_size_y=1;
 domain_size_z=12;

 mesh_x=0.01; 
 mesh_y=0.01;
 mesh_z=0.01;
 
 eps_zero=0.1;      % Value of homogenous media   
 sig_zero=0.009;      % value of homogenous media 
 E0 = 8.85418781762039080*1e-12;  % convert the abs eps to relative eps 
 

 
 time_window=200e-9;
 pml_number=20;
%

%% Calculate the discritization for New Starting models

nx_3=domain_size_x/mesh_x; 
nz_3=domain_size_z/mesh_z; 


%%  Read in true model

%__________ EPS __________% 
Fid1= fopen('model_syn_2016.eps','rb');
nz_1=fread(Fid1,1,'int');
nx_1=fread(Fid1,1,'int');
EV=fread(Fid1,nx_1*nz_1,'double');fclose(Fid1);
EV=reshape(EV,nx_1,nz_1);

%__________ SIG __________% 

Fid2= fopen('model_syn_2016.sig','rb');
nz_2=fread(Fid2,1,'int');
nx_2=fread(Fid2,1,'int');
SA=fread(Fid2,nx_2*nz_2,'double');fclose(Fid2);
SA=reshape(SA,nx_2,nz_2);


%__________ Pre-process data __________% 
                   
EV=EV/(8.85418781762039080*1e-12);


SA_1000=1000*SA;
SA_log=log10(SA_1000);


    
 %%  Creat the starting models 
 
eps_sm=sprintf('model_High_res.eps');
sig_sm=sprintf('model_High_res.sig');

Fid3= fopen (eps_sm, 'w');
Fid4= fopen (sig_sm, 'w');

fwrite(Fid3,nz_3,'int');
fwrite(Fid3,nx_3,'int');

fwrite(Fid4,nz_3,'int');
fwrite(Fid4,nx_3,'int');


eps_table=zeros(nx_2,nz_2); 
sig_table=zeros(nx_2,nz_2); 


%% Add layering or structure to starting models 


    for i=1:nx_3
       for j=nz_3

%_____________ Homogenous BG ______________%

eps_table(i,j)=eps_zero; 
sig_table(i,j)=sig_zero; 
                      
%____________ Layer & Object 1 ________________%            
%          if    i<67 && i>58 ;
%                 TM_Eps1(i,j)=20;
%          else    
%               TM_Eps1(i,j)=TM_Eps(i,j);
% 
%          end        
         
       end 
    end


  
  %% write down the EPS and SIG tables   
    
 %
 %    note : both need an converting units to be useful! 
 %  
  
 fwrite(Fid3,eps_table,'double');    
 fwrite(Fid4,sig_table,'double');    

fclose(Fid3);
fclose(Fid4); 

  
    
