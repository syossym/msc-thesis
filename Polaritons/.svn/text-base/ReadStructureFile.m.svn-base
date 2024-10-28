function structure = ReadStructureFile()

[filename, pathname, filterindex] = uigetfile({'*.str', 'Structure file'; '*.*',  'All Files'});
imported_file = importdata([pathname filename], '\t');
disp(['Imported structure file: ' filename]);

for (mm=1:length(imported_file.data))
   structure{mm}.Name = imported_file.textdata{6+mm,1};
   structure{mm}.x = imported_file.data(mm,1);
   structure{mm}.L = imported_file.data(mm,2);       % [A]
   structure{mm}.Active = imported_file.data(mm,3);   
   structure{mm}.N = imported_file.data(mm,4);       % [cm-^3]
end
