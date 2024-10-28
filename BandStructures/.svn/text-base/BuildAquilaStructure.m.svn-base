
% Read the structure file
Structure = ReadStructureFile();

% Build the AQUILA structure
width = 0;
for (ss=1:length(Structure))
    if (Structure{ss}.L > 100)
        res = 10;
    elseif (Structure{ss}.L < 100 && Structure{ss}.L > 50)
        res = 5;
    elseif (Structure{ss}.L < 50 && Structure{ss}.L > 10)
        res = 2;
    else
        res = 1;
    end
    
    if (Structure{ss}.Active)
        %disp(['Adding active layer: ' Structure{ss}.Name]);
        add_qbox([width, width+Structure{ss}.L], 1, 3, GE);
        %add_qbox([width, width+Structure{ss}.L], 1, 3, HH);
        add_pbox([width, width+Structure{ss}.L], CB+VB);
    end
    %disp(['Adding passive layer: ' Structure{ss}.Name]);
    add_mbox(Structure{ss}.L, 1, Structure{ss}.x, Structure{ss}.N);
    
    width = width + Structure{ss}.L;
end
add_pbox([0, width], CB+VB);
