function FID_mat = read_shape_bruker_data(path,Method_Params)

%Read data and shape in order of 
%[points x projections x Slices x echoes x bvalues x coils x repetitions]
try
    FIDs = Bruker_Load(fullfile(path,'rawdata.job0'));
catch
    try
        FIDs = Bruker_Load(fullfile(path,'fid'));
    catch
        FIDs = Bruker_Load(fullfile(path,'ser'));
    end
end

FID_mat = zeros(Method_Params.NPts,Method_Params.NPro,Method_Params.NSlices,Method_Params.NumTEs,Method_Params.Nbvalue,Method_Params.NCoil,Method_Params.Repetitions);

%I think the best way to do this will be to use the number of points to
%sort:
FIDs = reshape(FIDs,Method_Params.NPts,[]);

%Looping order:
%Coils inside of Slices, inside of bvalue, inside of echoes, inside of
%projections, inside of repetitions

looping = zeros(6,size(FIDs,2));

%This is probably horribly inefficient, but it's easy and relatively
%foolproof.
%looping: row 1 is coils
%         row 2 is slices
%         row 3 is N b-values
%         row 4 is N TEs
%         row 5 is Projections
%         row 6 is repetitions

count = 1;
for i = 1:Method_Params.Repetitions
    for j = 1:Method_Params.NPro
        for k = 1:Method_Params.NumTEs
            for l = 1:Method_Params.Nbvalue
                for m = 1:Method_Params.NSlices
                    for n = 1:Method_Params.NCoil
                        looping(1,count) = n;
                        looping(2,count) = m;
                        looping(3,count) = l;
                        looping(4,count) = k;
                        looping(5,count) = j;
                        looping(6,count) = i;
                        count = count+1;
                    end
                end
            end
        end
    end
end

for i = 1:size(FIDs,2)
    FID_mat(:,looping(5,i),looping(2,i),looping(4,i),looping(3,i),looping(1,i),looping(6,i)) = FIDs(:,i);
end




