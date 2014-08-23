function varargout = phxsim2(varargin)
% PHXSIM2 M-file for phxsim2.fig
%      PHXSIM2, by itself, creates a new PHXSIM2 or raises the existing
%      singleton*.
%
%      H = PHXSIM2 returns the handle to a new PHXSIM2 or the handle to
%      the existing singleton*.
%
%      PHXSIM2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHXSIM2.M with the given input arguments.
%
%      PHXSIM2('Property','Value',...) creates a new PHXSIM2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phxsim2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phxsim2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phxsim2

% Last Modified by GUIDE v2.5 08-Mar-2010 11:22:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @phxsim2_OpeningFcn, ...
    'gui_OutputFcn',  @phxsim2_OutputFcn, ...
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


% --- Executes just before phxsim2 is made visible.
function phxsim2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phxsim2 (see VARARGIN)

% Choose default command line output for phxsim2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes phxsim2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = phxsim2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

START = str2double(get(handles.inputStart,'String'));
END = str2double(get(handles.inputEnd,'String'));
Charge = str2double(get(handles.inputCharge,'String'));

table1 = get(handles.uitable1,'Data');
Temp=table1{1,2};
Time_p=table1{2,2}/1000;
pH_p=table1{3,2};
QFratio=[table1{4,2},table1{5,2},table1{6,2},table1{7,2}];

table2 = get(handles.uitable2,'Data');
pf_U=ones(1,size(table2,1));
pf_I=ones(1,size(table2,1));
pf_N=ones(1,size(table2,1));
clear Seq
for i=1:size(table2,1)
    pf_U(i)=table2{i,3};
    pf_I(i)=table2{i,4};
    pf_N(i)=table2{i,5};
    Seq(i)=table2{i,2};
end

table3 = get(handles.uitable3,'Data');
fractionU=table3{1,1}/100;
fractionI=table3{2,1}/100;
fractionN=table3{3,1}/100;
if fractionU+fractionI+fractionN~=1
    error('The sum of U, I & N fractions must be 100!')
end

% obsPeaks_U=phxsim2_sim(START,END,Charge,Temp,Time_p,pH_p,QFratio,pf_U,Seq); %call phxsim_sim.m
% obsPeaks_I=phxsim2_sim(START,END,Charge,Temp,Time_p,pH_p,QFratio,pf_I,Seq);
% obsPeaks_N=phxsim2_sim(START,END,Charge,Temp,Time_p,pH_p,QFratio,pf_N,Seq);
%%%2010-04-27 changed to calling hxsim.m:
fractionD=QFratio(1)/(QFratio(1)+QFratio(2)+QFratio(3));
obsPeaks_U=hxsim(START,END,Charge,Seq,pf_U,Temp,pH_p,fractionD,Time_p,1);
obsPeaks_I=hxsim(START,END,Charge,Seq,pf_I,Temp,pH_p,fractionD,Time_p,1);
obsPeaks_N=hxsim(START,END,Charge,Seq,pf_N,Temp,pH_p,fractionD,Time_p,1);

clear obsPeaks
obsPeaks(:,1)=obsPeaks_U(:,1);
obsPeaks(:,2)=fractionU*obsPeaks_U(:,2) + fractionI*obsPeaks_I(:,2) + fractionN * obsPeaks_N(:,2);

[fractionFitU, fractionFitN, fractionFitI, DistrFitI]=phx_cfit(obsPeaks_U(:,2), obsPeaks_N(:,2), obsPeaks(:,2)); %call phx_cfit.m
            
set(handles.axes1,'NextPlot','replace')
stem(handles.axes1,obsPeaks(:,1),obsPeaks(:,2),'k')

set(handles.axes1,'NextPlot','add')
plot(handles.axes1,obsPeaks(:,1), fractionU*obsPeaks_U(:,2),'m:');
set(handles.axes1,'NextPlot','add')
plot(handles.axes1,obsPeaks(:,1), fractionN*obsPeaks_N(:,2),'b:');
set(handles.axes1,'NextPlot','add')
plot(handles.axes1,obsPeaks(:,1), fractionI*obsPeaks_I(:,2),'g:');

set(handles.axes1,'NextPlot','add')
plot(handles.axes1,obsPeaks(:,1), fractionFitU*obsPeaks_U(:,2),'m','LineWidth',1.5);
set(handles.axes1,'NextPlot','add')
plot(handles.axes1,obsPeaks(:,1), fractionFitN*obsPeaks_N(:,2),'b','LineWidth',1.5);
set(handles.axes1,'NextPlot','add')
plot(handles.axes1,obsPeaks(:,1), DistrFitI,'g','LineWidth',1.5);

table3{1,2}=fractionFitU*100;
table3{2,2}=fractionFitI*100;
table3{3,2}=fractionFitN*100;
set(handles.uitable3,'Data',table3);

% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function inputStart_Callback(hObject, eventdata, handles)
% hObject    handle to inputStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputStart as text
%        str2double(get(hObject,'String')) returns contents of inputStart as a double


% --- Executes during object creation, after setting all properties.
function inputStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputEnd_Callback(hObject, eventdata, handles)
% hObject    handle to inputEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputEnd as text
%        str2double(get(hObject,'String')) returns contents of inputEnd as a double


% --- Executes during object creation, after setting all properties.
function inputEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputCharge_Callback(hObject, eventdata, handles)
% hObject    handle to inputCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputCharge as text
%        str2double(get(hObject,'String')) returns contents of inputCharge as a double


% --- Executes during object creation, after setting all properties.
function inputCharge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)





