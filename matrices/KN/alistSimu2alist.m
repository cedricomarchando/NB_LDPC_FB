% alist2matrix_nb.m
% author: Cédric Marchand
% construct parity check matrix then build alist matrix

clear all;
clc;




fileID =fopen('N5760_K1920_GF64_ali.txt');


N = fscanf(fileID,'%d',1);
M = fscanf(fileID,'%d',1);
GF = fscanf(fileID,'%d',1);
%dvmax = fscanf(fileID,'%d',1);
%dcmax = fscanf(fileID,'%d',1);

dc=zeros(1,M);
dv=zeros(1,N);

matrix=zeros(M,N);
matrix=matrix-1;
matrix2=matrix;

for i=1:N
   dv(i) = fscanf(fileID,'%d',1); 
end

for i=1:M
   dc(i) = fscanf(fileID,'%d',1); 
end


for i=1:M
    for j=1:dc(i)
       pos = fscanf(fileID,'%d',1); 
       value = fscanf(fileID,'%d',1);
       matrix2(i,pos)=value;
    end
end


fclose(fileID);






H = matrix2;

H_bin = H;
H_bin( H > -1 ) = 1;
H_bin( H<0 ) = 0;


[M,N] = size(H);
dc = sum(H_bin,2);
dv = sum(H_bin,1);
dcmax = max(dc);
dvmax = max(dv);

[nlist,~] = find(H_bin);
[mlist,~] = find(H_bin');

GF=64;

%write in a file named alist.txt
     fileID =fopen('alist.txt','at');    
      fprintf(fileID,'%d %d %d \n',N,M,GF);
    for i=1:N
        fprintf(fileID,'%d ',dv(i));
    end
    fprintf(fileID,' \n');
    for i=1:M
        fprintf(fileID,'%d ',dc(i));
    end
    fprintf(fileID,' \n');
    k=0;
    for i=1:N
        for j=1:dv(i)
            fprintf(fileID,'%d ',nlist(k + j ));
            fprintf(fileID,'%d ',H( nlist(k + j  ),i));
        end
        k=k+dv(i);
        fprintf(fileID,' \n');
    end
    
    fprintf(fileID,' \n');
    
    k=0;
    for i=1:M
        for j=1:dc(i)
            fprintf(fileID,'%d ',mlist(k + j ));
            fprintf(fileID,'%d ',H( i, mlist(k + j  )));
        end
        k=k+dc(i);
        fprintf(fileID,' \n');
    end

    
      fclose(fileID);   

%end



