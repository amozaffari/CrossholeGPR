%clc; clear ;

% ____________ Reading ________________ %%
    
fid=fopen('source.dat','rb');% fopen ('source.dat','r','b') um int32 zulesen
length=fread (fid,1,'int');
Src=fread(fid,(length),'double');
dt=fread(fid,1,'double');
fclose (fid);





%%%_______________ Writing _______________ %%


% Fid1= fopen('source_resampled741.dat','wb');
% 
% 
% fwrite(Fid1,length,'int');  % number of sampels 
% 
% for sample_num=1:length
%  fwrite(Fid1,Src(sample_num),'double');                                               
% end
% fwrite(Fid1,dt ,'double');  % number of sampels    
% 
%                                             
% fclose(Fid1);
                                            
