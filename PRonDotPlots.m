% Author: Konstantinos Prousalis
% Email: kprousalis@csd.auth.gr

% PURPOSE
% This program is designed to detect imperfect diagonal lines for dot matrix
% plots. As accepted files are only fasta format.
% Only pairwise alignment is allowed for different sequences covering the 
% most common patterns:
% 1. perfect match
% 2. homologous in the same diagonal
% 3. homologous that are not in the same diagonal
% 4. indels by requesting via dotplot
% Plotting a sequence by itself may form patterns that are not supported.

% HINT 1
% The user can download the readily available fasta files from genomic 
% databases, save them as querySeq.fasta and refSeq.fasta and run their content.
% Beware, remove some initial instructions.

% HINT 2
% The user can use the seqdotplot(seq1,seq2) function throughout the code
% to enable the view of the dotplot per window.

% HINT 3
% Runtime may suffer an exponential slow down for very large windows 
% (better to use windows smaller than 1000) or wide scanning around the 
% main diagonal (critical parameter hw, see the code).


function varargout = PRonDotPlots(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PRonDotPlots_OpeningFcn, ...
                   'gui_OutputFcn',  @PRonDotPlots_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before PRonDotPlots is made visible.
function PRonDotPlots_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

dataSeq = struct('querySeq',0,'refSeq',0);
set(handles.pushbutton1,'UserData',dataSeq);
dataWin = struct('windowsize',0);
set(handles.popupmenu1,'UserData',dataWin);
dataSence = struct('sensitivity',0);
set(handles.popupmenu2,'UserData',dataSence);


% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = PRonDotPlots_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)

str=get(hObject,'String');
val=get(hObject,'Value');
wsize=str2double(str{val});

dataWin = get(hObject,'UserData');
dataWin.windowsize=wsize;
set(hObject,'UserData',dataWin);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

global totalMtrx;
global maxQueryPieces;
global maxRefPieces;

totalMtrx=zeros;

%remove newline delimeter from texts
filecontent1 = fileread('refSeq.fasta');
newcontent1 = regexprep(filecontent1,'\s+','');
fid = fopen('refSeq.fasta', 'w');
fwrite(fid, newcontent1);
refSeqSize = size(filecontent1);
fclose(fid);

filecontent2 = fileread('querySeq.fasta');
newcontent2 = regexprep(filecontent2,'\s+','');
fid = fopen('querySeq.fasta', 'w');
fwrite(fid, newcontent2);
querySeqSize=size(filecontent2);
fclose(fid);

dataWin = get(handles.popupmenu1,'UserData');
%dataWin.windowsize
dataSense = get(handles.popupmenu2,'UserData');
%dataSense.sensitivity

maxQueryPieces=floor(querySeqSize(2)/dataWin.windowsize);
maxRefPieces=floor(refSeqSize(2)/dataWin.windowsize);
width=floor((dataSense.sensitivity/100)*dataWin.windowsize);
threshold75=0.75*dataWin.windowsize;
threshold50=0.50*dataWin.windowsize;
threshold25=0.25*dataWin.windowsize;

counter75=0;
counter50=0;
counter25=0;

fid1 = fopen('querySeq.fasta', 'r');
fid2 = fopen('refSeq.fasta', 'r');

% CRITICAL PARAMETER 1
% Set the scanning width along the main diagonal in terms of windows
% this parameter will slow down the execution speed
hw=0;
mainDiagWidth=2*hw+1;

% CRITICAL PARAMETER 2
% Set the starting point in the frame of the overall dotplot
offsetQwr=0;
offsetRef=0;
fseek(fid1,offsetQwr,-1);
fseek(fid2,offsetRef,-1);


for c = 1:175%maxQueryPieces
    seqQwr=fscanf(fid1,'%c',dataWin.windowsize);
    %frewind(fid2);
    %fseek(fid2,offsetRef,-1);
    if c < hw
        %frewind(fid2);
        fseek(fid2,offsetRef,-1);
    else
        fseek(fid2,((c-hw)*dataWin.windowsize),-1);
    end
    for r = 1:mainDiagWidth%maxRefPieces
        seqRef=fscanf(fid2,'%c',dataWin.windowsize);
        [Matches, Matrix] = seqdotplot(seqQwr, seqRef);close(1);
        mtrx = full(Matrix);
        for i = -(width):(width)
            d = nnz(diag(mtrx,i));
            if d > threshold75
                counter75=counter75+1;
                X = ['Read: @ref x = ',num2str(c),' and @query y = ',num2str(r),' is positive 75%'];disp(X);
                if c-1 < hw
                    totalMtrx(c,r) = 75;
                else
                    totalMtrx(c,c-hw+r) = 75;
                end
            elseif d > threshold50
                counter50=counter50+1;
                X = ['Read: @ref x = ',num2str(c),' and @query y = ',num2str(r),' is positive 50%'];disp(X);
                if c-1 < hw
                    totalMtrx(c,r) = 50;
                else
                    totalMtrx(c,c-hw+r) = 50;
                end
            end
        end
    end
end
fclose(fid2);
fclose(fid1);

% Result information
X1 = ['Query  seq    size: ',num2str(querySeqSize(2))];
X2 = ['Reference seq size: ',num2str(refSeqSize(2))];
X3 = ['Window        size: ',num2str(dataWin.windowsize)];
X4 = ['Window    coverage: ',num2str(dataSense.sensitivity),'% or width ',num2str(width),' number of diagonals up or down.'];
X5 = ['Total Query pieces: ',num2str(maxQueryPieces)];
X6 = ['Total Ref   pieces: ',num2str(maxRefPieces)];
disp(X1);disp(X2);disp(X3);disp(X4);disp(X5);disp(X6);
counter=counter75+counter50;
R = ['Cutoff 50%: ',num2str(counter50),' and cutoff 75%: ',num2str(counter75),'         Total hits: ',num2str(counter)];
disp(R);
disp("Scanning is finished!");




% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)

str=get(hObject,'String');
val=get(hObject,'Value');
sensitivity=str2double(str{val});

dataSense = get(hObject,'UserData');
dataSense.sensitivity=sensitivity;
set(hObject,'UserData',dataSense);

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

global totalMtrx;
global maxQueryPieces;
global maxRefPieces;


figure
% May be more convinient to rotate or flip your results
% example... imagesc(rot90(rot90(rot90(flipud(totalMtrx)))))
imagesc(totalMtrx)
colormap(flipud(gray))
colorbar
colorbar('YTick',[0 10 20 30 40 50 60 70 80 90 100])
caxis([0 100])
axis square
xticks('auto')
yticks('auto')
title('Map of positive windows')
xlabel('Reference Pieces')
ylabel('Query Pieces')