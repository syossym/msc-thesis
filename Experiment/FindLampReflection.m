function LampRef = FindLampReflection(LampRefs, offset)

LampRef = [];

for (ii=1:length(LampRefs))
   if (strcmp(LampRefs{ii}.Spectral_Offset, offset))
       LampRef = LampRefs{ii}.Data(:,2);
   end
end

if (isempty(LampRef))
    LampRef = ones(size(LampRefs{1}.Data(:,2)));
end
