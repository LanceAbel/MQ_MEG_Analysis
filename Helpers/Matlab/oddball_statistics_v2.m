
windows_machine = ispc
if windows_machine
    details_file = 'E:\subject.txt';
else
    details_file = '/subject.txt'
end
fid=fopen(details_file);
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);
subject = char(tlines(1));
type = char(tlines(2));
data_root = char(tlines(4));
data_root_mac = char(tlines(5));
matlab_general_code = char(tlines(6));
code_folder = char(tlines(7));
code_folder_mac = char(tlines(8));
NUM_TONES =  str2double(tlines(16));
tones = 1:NUM_TONES
tone_repetitions = 1:NUM_TONES;
% If running on single person, uncomment
if ~windows_machine
    code_folder = code_folder_mac
end


cd(matlab_general_code)
subjects = [];
folder_list = workspace_find_variable_list_files(strcat(data_root,'\'));
%disp(folder_list);
num_subjects = 0
for i=1:length(folder_list)
    folder = char(folder_list(i));
    regex = regexp(folder, '(\d\d\d\d)', 'match'); % Only folders containing subject codes only
    try
        subj = (regex(1));
        subjects = [subjects, subj];
    end
end
subjects = unique(subjects);
disp(subjects)

rm_subjects = {'3434','2607'}
for j=1:length(rm_subjects)
    subjects(ismember(subjects,rm_subjects(j))) = [];
end




% % Import the file
% newData1 = load('-mat', fileToRead1);
% 
% % Create new variables in the base workspace from those fields.
% vars = fieldnames(newData1);
% for i = 1:length(vars)
%     assignin('base', vars{i}, newData1.(vars{i}));
% end




% % Counts

tone_counts = struct('count', num2cell(tones))
tone_repetition_counts = struct('count', num2cell(tone_repetitions))

tone_counts_all = []

% Calculate the statistics of tones played
for j =1:length(subjects)

    folder = strcat(data_root,subjects(j),'\'); % folder with subject subfolders
    folder = char(folder);
    cd(folder)

    disp(subjects(j))

    if type == 'Child'
        all_txt = dir('*oddball*.txt');
        num_text_files = length(all_txt);
        if num_text_files > 0
            all_txt_table = struct2table(all_txt); % convert the struct array to a table
            all_txt_by_size = sortrows(all_txt_table, 'bytes', 'descend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
            all_txt_by_size = table2struct(all_txt_by_size);
            txtfile = all_txt_by_size;    % Gets first one somehow
            txtfile = txtfile(1);
            txtfile_name = fullfile(txtfile.folder,txtfile.name);
        
            if contains(txtfile_name, 'oddball')
                %disp(txtfile_name
                %fid = fopen(txtfile_name);
                data = load(txtfile_name);
                tone_data = data(:,3);
  
                tone_repetitions = 0;
                for L = 1:length(tone_data)
                    tone_val = tone_data(L);
                    if L > 1
                        tone_val_prev = tone_data(L-1);
                        if tone_val_prev == tone_val
                            tone_repetitions = tone_repetitions + 1;
                        else
                            tone_repetition_counts(tone_repetitions).count = tone_repetition_counts(tone_repetitions).count + 1;
                            tone_repetitions = 1;
                        end
                    end
                    %Assign value
                    tone_counts(tone_val).count = tone_counts(tone_val).count + 1;
                end
   
            end
        end
        

     elseif type == 'Adult'
        try
            all_txt = dir('*event_meg.txt');
            num_text_files = length(all_txt);
            if num_text_files > 0
                all_txt_table = struct2table(all_txt); % convert the struct array to a table
                all_txt_by_size = sortrows(all_txt_table, 'bytes', 'descend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
                all_txt_by_size = table2struct(all_txt_by_size);
                txtfile = all_txt_by_size;    % Gets first one somehow
                txtfile = txtfile(1);
                txtfile_name = fullfile(txtfile.folder,txtfile.name);
    
    
                disp(subjects(j))
                % Skip first line
                [pathstr, name, ext] = fileparts(txtfile_name);         %had to add this line to get the file extension
                if strcmp(ext, 'txt')
                   data = importdata(txtfile_name, ',', 1);              %I am guessing about the delimiter
                else
                   data = importdata(txtfile_name);
                end
                data = data.data;
                %data = load(txtfile_name);
    
                tone_data = data(:,2);
    
                tone_repetitions = 0;
                for L = 1:length(tone_data)
                    tone_val = tone_data(L);
                    if L > 1
                        tone_val_prev = tone_data(L-1);
                        if tone_val_prev == tone_val
                            tone_repetitions = tone_repetitions + 1;
                        else
                            tone_repetition_counts(tone_repetitions).count = tone_repetition_counts(tone_repetitions).count + 1;
                            tone_repetitions = 1;
                        end
                    end
                    %Assign value
                    tone_counts(tone_val).count = tone_counts(tone_val).count + 1;
                end
            end
        catch
            disp("Issue ")
        end

    end
end



%1620/18157

sums_ = 0;
for j = 1:length(tone_counts)
    sums_ = sums_ + tone_counts(j).count;
end
disp(sums_)


sums_ = 0;
for j = 1:length(tone_repetition_counts)
    sums_ = sums_ + tone_repetition_counts(j).count;
end
disp(sums_)





%% UNFORTUNATELY THE DISTRIBUTION OF REPETITIONS IS DIFFERENT FOR CHILDREN AND ADULTS

