function writeL2C(filename,L2valdata)
%%
fid = fopen([filename 'fAPAR.csv'],'w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,\r\n','time','sza','fAPAR','fAPAR_unc','fAPARchl','fAPARchl_unc');
for k = 1:length(L2valdata)
    time = (L2valdata(k).results.L2biophys.t);
    for j = 1:length(L2valdata(k).results.L2biophys.fAPAR)
        %
        x = [L2valdata(k).results.L2biophys.sza(j), ...
            L2valdata(k).results.L2biophys.fAPAR(j), ...
            L2valdata(k).results.L2biophys.fAPAR_unc, ...
            L2valdata(k).results.L2biophys.fAPARchl(j), ...
            L2valdata(k).results.L2biophys.fAPARchl_unc];

        %fprintf(fid,'%20s,%5.4f,%6.4f,%6.4f,%6.4,%6.4f\r\n',time(j,:),x);
        fprintf(fid,'%20s',time(j,:));
        for z = 1:length(x)
            fprintf(fid,',%5.4f',x(z));
        end
         fprintf(fid,'\r\n');

    end
end
fclose(fid);
%%
fid = fopen([filename 'APARchl.csv'],'w');
fprintf(fid,'%s,%s,%s\r\n','time','APARchl','APARchl_unc');
for k = 1:length(L2valdata)
    time = (L2valdata(k).results.L2biophys.time_full);
    for j = 1:length(L2valdata(k).results.L2biophys.APARchl)
        fprintf(fid,'%20s,%7.2f,%7.2f\n',time(j,:),L2valdata(k).results.L2biophys.APARchl(j),L2valdata.results.L2biophys.APARchl_unc(j));
    end
end
fclose(fid);
%%
fid = fopen([filename 'bio.csv'],'w');
   fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\r\n','time','LCC','LCC_unc','LCCAR','LCCAR_unc','LAR','LAI_unc');

for k = 1:length(L2valdata)
time = (L2valdata(k).results.L2biophys.time_full);

        fprintf(fid,'%20s,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f\n',time(1,1:11),L2valdata(k).results.L2biophys.LCC,L2valdata.results.L2biophys.LCCunc, ...
            L2valdata(k).results.L2biophys.LCAR,L2valdata.results.L2biophys.LCARunc,...
            L2valdata(k).results.L2biophys.LAI,L2valdata.results.L2biophys.LAIunc);
end
fclose(fid);

%%
fid = fopen([filename 'fesc.csv'],'w');
fprintf(fid,'%s,%s,%s\r\n','time','sza','fesc(spectrum,642:1:846nm)');
for k = 1:length(L2valdata)
    time = (L2valdata(k).results.L2biophys.t);
    for j = 1:size(time,1)
        fprintf(fid,'%20s,%5.4f',time(j,:),L2valdata(k).results.L2biophys.sza(j));

        for z = 1:205
            fprintf(fid,',%5.4f',L2valdata(k).results.L2biophys.sigmaF(z+2,j));
        end
        fprintf(fid,'\r\n');
    end
end
fclose(fid);

%%
fid = fopen([filename 'fesc_unc.csv'],'w');
fprintf(fid,'%s,%s,%s\r\n','time','fesc(spectrum,642:1:846nm)');
for k = 1:length(L2valdata)
    time = (L2valdata(k).results.L2biophys.t);
    fprintf(fid,'%20s',time(1,:));
    for z = 1:205
        fprintf(fid,',%5.4f',L2valdata(k).results.L2biophys.sigmaF_unc(z+2));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);


if ~isfield(L2valdata(k).results.L2biophys, 'FQE')
fid = fopen([filename 'FQE.csv'],'w');
fprintf(fid,'%s,%s,%s\r\n','time','FQE','FQE_unc');
for k = 1:length(L2valdata)
    time = (L2valdata(k).results.L2biophys.time_full);
    for j = 1:length(L2valdata(k).results.L2biophys.APARchl)
        fprintf(fid,'%20s,%7.2f,%7.2f\n',time(j,:),L2valdata(k).results.L2biophys.FQE(j),L2valdata.results.L2biophys.FQE_unc(j));
    end
end
fclose(fid);
end
