clear all
EV=[1 3 5;7 10 6;1 2 3];

num_inter=2; %number of interpolations
num_inter_1=num_inter+1; %number of interpolations +1

EV_interpolated=zeros(size(EV,1)+num_inter*(size(EV,1)-1),size(EV,2)+num_inter*(size(EV,2)-1));

for i=1:size(EV,1)
    for j=1:size(EV,2)
        EV_interpolated(num_inter_1*(i-1)+1,num_inter_1*(j-1)+1)=EV(i,j);
    end
end

%%

% row interpolation

for i=1:size(EV,1)
    for j=1:size(EV,2)-1
        for k=1:num_inter
            EV_interpolated(num_inter_1*(i-1)+1,num_inter_1*(j-1)+1+k)=EV(i,j)+k*(EV(i,j+1)-EV(i,j))/num_inter_1;
        end
    end
end

%%

% column interpolation

for j=1:size(EV,2)
    for i=1:size(EV,1)-1
        for k=1:num_inter
            EV_interpolated(num_inter_1*(i-1)+1+k,num_inter_1*(j-1)+1)=EV(i,j)+k*(EV(i+1,j)-EV(i,j))/num_inter_1;
        end
    end
end

%%

% fill in the rest

for i=1:size(EV_interpolated,1)
    for j=1:size(EV,2)-1
        for k=1:num_inter
            EV_interpolated(i,num_inter_1*(j-1)+1+k)=EV_interpolated(i,num_inter_1*(j-1)+1)+k*(EV_interpolated(i,j*num_inter_1+1)-EV_interpolated(i,num_inter_1*(j-1)+1))/num_inter_1;
        end
    end
end


%%
