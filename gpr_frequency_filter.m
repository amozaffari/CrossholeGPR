function [gprdat, para]=gpr_frequency_filter(indata, para)

gprdat=indata;



try  % check if ntaper defined
              para.freq.ntaper= para.freq.ntaper;
catch
          para.freq.ntaper = 11;
 end


for fn=1:size(indata,2)

    fdata=indata(1,fn); 
    
[nsamp ntrace]=size(fdata.data);
nf=2^ceil(log2(fdata.no_samples));
fdatag=zeros(nf,fdata.no_traces);
fdatag(1:fdata.no_samples,:)=fdata.data;
fdatag=fft(fdatag,nf);
df=1/(nf*fdata.dt*1e-9);  %step freq
f=(0:(nf-1))*df;            % x axis

mfdatag  = mean(abs(fdatag(:,min(~isnan(fdatag)))),2);  % for plotting takes mean

if(~isfield(para.freq,'corner_frequencies')) 
    figure;
    hold on
    plot(f(1:end/2)*1e-6,mfdatag(1:end/2));
    hold on 
    plot([fdata.freq fdata.freq],[0 max(mfdatag(1:end/2))], 'black')
   xlim([0 3*fdata.freq])
    bf = inputdlg('Which is the value of the starting frequency?','',1);
    fs = str2num(bf{:});
    ef = inputdlg('Which is the value of the ending frequency?','',1);
    fe = str2num(ef{:});
    close 
else
    fs=para.freq.corner_frequencies(1);
    fe=para.freq.corner_frequencies(2);
end  
    bfs=round(fs*1e6/df)+1;
    efs=round(fe*1e6/df)+1;
%ntap2=round((para.freq.ntaper-1)/2);  
ntap2=((para.freq.ntaper-1)/2); % check if filter works
if bfs<ntap2-1
    error(['Taper is too long, should be less then ' num2str(2*bfs-1)])
end
        
    taper='cosine';
if taper=='cosine' ;
    tap=(pi).*(0:(para.freq.ntaper-1))/(para.freq.ntaper-1);
    itap=transpose((cos(tap).'+1)/2);
    tap=fliplr(itap);
elseif taper=='linear';
    ntap2=(para.freq.ntaper-1)/2;
    tap=(0:para.freq.ntaper-1)/8;
    itap=fliplr(tap);
end
fdatagf=fdatag;

for itrace=1:fdata.no_traces
    fdatagf(bfs-ntap2:bfs+ntap2,itrace)=transpose(tap).*fdatag(bfs-ntap2:bfs+ntap2,itrace);
    fdatagf(nf-bfs-ntap2:nf-bfs+ntap2,itrace)=transpose(itap).*fdatag(nf-bfs-ntap2:nf-bfs+ntap2,itrace);
end
fdatagf(1:bfs-ntap2,:)=0;
fdatagf(efs+ntap2:nf-efs-ntap2,:)=0;
fdatagf(nf-bfs+ntap2:nf,:)=0;
for itrace=1:fdata.no_traces
    fdatagf(efs-ntap2:efs+ntap2,itrace)=transpose(itap).*fdatag(efs-ntap2:efs+ntap2,itrace);
    fdatagf(nf-efs-ntap2:nf-efs+ntap2,itrace)=transpose(tap).*fdatag(nf-efs-ntap2:nf-efs+ntap2,itrace);
end
idatagf=real(ifft(fdatagf,nf));


gprdat(1,fn).data = idatagf(1:fdata.no_samples,:);


    
    mfdatagf = mean(abs(fdatagf(:,min(~isnan(fdatag)))),2);

   fig= figure;
    hold on
    plot(f(1:end/2), mfdatag(1:end/2));
    plot(f(1:end/2), mfdatagf(1:end/2),'r');
    xlabel('Frequency (GHz)')
    legend('complete', 'used')
    
      
   name_tmp1 = [para.output_path gprdat(1,fn).name 'freq_spectrum'  para.plt.fmt];    
   saveas (fig,name_tmp1,'png') 
   close all
                     


end 
para.freq.corner_frequencies = [fs fe];
      %% plot  data 
 gpr_plot_1fg( gprdat,para, 'freq-filter' )

 
end