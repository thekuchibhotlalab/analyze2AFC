% Create a new figure window
fig = figure('Position',[100 100 800 600],'Name','Data Analysis GUI');

% Add a button to open a file selection dialog box
uicontrol('Parent',fig,'Style','pushbutton','String','Open Files',...
          'Position',[20 20 150 30],'Callback',@openFiles);

% Add a selection menu for different functions
uicontrol('Parent',fig,'Style','popupmenu','String',{'Mean','Median','Mode'},...
          'Position',[180 20 150 30],'Callback',@applyFunction);

% Add a text area to display the data and results
uicontrol('Parent',fig,'Style','text','Position',[20 60 760 520],'String','',...
          'HorizontalAlignment','left','FontSize',12);

% Callback function for the "Open Files" button
function openFiles(src,event)
    % Open a file selection dialog box and load the selected files
    [file1,path1] = uigetfile('*.csv');
    [file2,path2] = uigetfile('*.txt');
    data1 = readtable(fullfile(path1,file1));
    data2 = readfile(fullfile(path2,file2));
    
    % Display the data from the files in the text area
    set(src,'String',{'Data from CSV table:',data1,...
                      'Data from text file:',data2});
end

% Callback function for the selection menu
function applyFunction(src,event)
    % Get the selected function from the menu
    func = get(src,'String');
    func = func{get(src,'Value')};
    
    % Apply the selected function to the data and display the result in the
    % text area
    switch func
        case 'Mean'
            % Calculate the mean of the data and display the result
            result = mean(data1);
            set(src,'String',{'Data from CSV table:',data1,...
                              'Data from text file:',data2});

    end
end 
