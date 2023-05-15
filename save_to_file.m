function [] = save_to_file(location,mean_aRMSE,aRMSE_odch,mean_Jx,...
    Jx_odch,mean_Jy,Jy_odch,mean_eps,eps_odch,cycles)

fid = fopen(append(location,'/aRMSE.csv'),'at');
fprintf(fid,'%f ',mean_aRMSE);
fprintf(fid,'%f ',aRMSE_odch);
fprintf(fid,'%f ',cycles);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(append(location,'/Jx.csv'),'at');
fprintf(fid,'%f ',mean_Jx);
fprintf(fid,'%f ',Jx_odch);
fprintf(fid,'%f ',cycles);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(append(location,'/Jy.csv'),'at');
fprintf(fid,'%f ',mean_Jy);
fprintf(fid,'%f ',Jy_odch);
fprintf(fid,'%f ',cycles);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(append(location,'/eps.csv'),'at');
fprintf(fid,'%f ',mean_eps);
fprintf(fid,'%f ',eps_odch);
fprintf(fid,'%f ',cycles);
fprintf(fid,'\n');
fclose(fid);

end