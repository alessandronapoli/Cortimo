
%-----------------Version 1.0 -----------------------%

%Dr Obeid gave me this function when I started acquiring EEG signals
%using the Modular Machine (OPENEEG). It loads EDF files by creating
%a specific structure. The data are saved in the matrix fd.data, this
%matrix must be reashaped in order to achieve useful data channels.


function fd = read_edf_file(fname)
%fd = file data

fid = fopen(fname,'rb');
fd.ver = fscanf(fid,'%c',8);
fd.loc_pt_info = fscanf(fid,'%c',80);
fd.loc_rc_info = fscanf(fid,'%c',80);
fd.date_start = fscanf(fid,'%c',8);
fd.time_start = fscanf(fid,'%c',8);
fd.nBytes_header = fscanf(fid,'%c',8);
fd.reserved1 = fscanf(fid,'%c',44);
fd.nDataRec = str2double(fscanf(fid,'%c',8));
fd.duration = str2double(fscanf(fid,'%c',8));
fd.ns = str2double(fscanf(fid,'%c',4));
ns = fd.ns;
fd.labels = fscanf(fid,'%c',[16,ns])';
fd.transducer = fscanf(fid,'%c',[80,ns])';
fd.physDim = fscanf(fid,'%c',[8,ns])';
fd.physMin = fscanf(fid,'%c',[8,ns])';
fd.physMax = fscanf(fid,'%c',[8,ns])';
fd.digiMin = fscanf(fid,'%c',[8,ns])';
fd.digiMax = fscanf(fid,'%c',[8,ns])';
fd.preFilter = fscanf(fid,'%c',[80,ns])';
fd.nSamplesPerRec = str2num(fscanf(fid,'%c',[8,ns])');
fd.reserved2 = fscanf(fid,'%c',[32,ns])';
fd.data = fread(fid,'short');
fclose(fid);
