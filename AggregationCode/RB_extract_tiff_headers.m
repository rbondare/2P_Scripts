function RB_extract_tiff_headers(recordingList, overwrite_intermediate)

dataDir   = 'Z:\Group Members\Rima\TEST\';
interBase = 'C:\Users\rbondare\IntermediateDir\';

for r = 1:numel(recordingList)
    recordingName = recordingList{r};
    basedir = fullfile(dataDir, recordingName);
    
    nameParts = strsplit(recordingName, '_');
    animalName = nameParts{1};
    
    if ~exist(basedir, 'dir')
        fprintf('Recording folder not found: %s\n', recordingName);
        continue;
    end
    
    interDir = fullfile(interBase, animalName);
    if ~exist(interDir,'dir'); mkdir(interDir); end
    IntermediateFilename = fullfile(interDir, [recordingName '_S2Presult.mat']);
    
    files_tif = dir(fullfile(basedir, '*.tif'));
    files_hdr = dir(fullfile(basedir, '*_header.mat'));
    if ~isempty(files_tif)
        files = files_tif;
    else
        files = files_hdr;
    end
    numFiles = numel(files);
    
    if numFiles == 0
        fprintf('No files found in %s\n', recordingName);
        continue;
    end
    
    Data = matfile(IntermediateFilename,'Writable',true);
    
    if ~exist(IntermediateFilename,'file') || overwrite_intermediate
        Headers = struct();
        w = waitbar(0, sprintf('Extracting headers: %s', recordingName));
        
        for f = 1:numFiles
            fullname = fullfile(files(f).folder, files(f).name);
            [~, file_name, ext] = fileparts(fullname);
            
            if strcmp(ext,'.tif')
                headerFile = fullfile(basedir, [file_name '_header.mat']);
                if ~exist(headerFile,'file')
                    [header,~,ImgInfo] = fast_get_scanimage_header_notMROI(fullname,'slice',1,'volume',1);
                    if ~isempty(header)
                        save(headerFile,'header','ImgInfo','-v7.3');
                    end
                else
                    load(headerFile,'header','ImgInfo');
                end
            else
                load(fullname,'header','ImgInfo');
            end
            
            Headers(f).ImgInfo = ImgInfo;
            if ~isempty(header)
                FN = fieldnames(header);
                for n = 1:numel(FN)
                    Headers(f).(FN{n}) = header.(FN{n});
                end
            end
            waitbar(f/numFiles, w)
        end
        close(w)
        
        k1 = find(~cellfun(@isempty,{Headers.SI}),1,'first');
        if ~isempty(k1)
            SI_Info = Headers(k1).SI;
            SI_Info.Clockstart = datetime(Headers(k1).epoch{1});
            
            if ~isfield(Headers(end).ImgInfo,'numVolumes') && Headers(end).SI.hStackManager.numSlices > 1
                numS = Headers(end-1).ImgInfo.numSlices + SI_Info. hFastZ.numDiscardFlybackFrames;
                FN = fieldnames(Headers);
                last = numel(Headers(end).frameNumbers) - rem(numel(Headers(end).frameNumbers), numS);
                for n = 3:numel(FN)
                    try
                        Headers(end).(FN{n}) = Headers(end).(FN{n})(1:last);
                    catch
                    end
                end
            end
            
            Data.Headers = Headers;
            Data.SI_Info = SI_Info;
            fprintf('Headers extracted for %s\n', recordingName);
        end
    end
end
end