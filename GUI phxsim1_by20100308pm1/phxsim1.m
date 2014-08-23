function varargout = phxsim1(varargin)
% PHXSIM1 M-file for phxsim1.fig
%      PHXSIM1, by itself, creates a new PHXSIM1 or raises the existing
%      singleton*.
%
%      H = PHXSIM1 returns the handle to a new PHXSIM1 or the handle to
%      the existing singleton*.
%
%      PHXSIM1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHXSIM1.M with the given input arguments.
%
%      PHXSIM1('Property','Value',...) creates a new PHXSIM1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phxsim1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phxsim1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phxsim1

% Last Modified by GUIDE v2.5 08-Mar-2010 10:19:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @phxsim1_OpeningFcn, ...
    'gui_OutputFcn',  @phxsim1_OutputFcn, ...
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


% --- Executes just before phxsim1 is made visible.
function phxsim1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phxsim1 (see VARARGIN)

% Choose default command line output for phxsim1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes phxsim1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = phxsim1_OutputFcn(hObject, eventdata, handles)
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
pf=ones(1,size(table2,1));
clear Seq
for i=1:size(table2,1)
    pf(i)=table2{i,4};
     Seq(i)=table2{i,2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%below part can be written as a separate function:

%%%data prep:
[peptideMass, distND, maxND, maxD]=pepinfo(Seq(START:END)); %call pepinfo.m
kcDH = fbmme_dh(Seq, pH_p, Temp, 'poly'); %call fbmme_dh.m
kcHD = fbmme_hd(Seq, pH_p, Temp, 'poly', 0); %call fbmme_hd.m
Dfraction=QFratio(1)/(QFratio(1)+QFratio(2)+QFratio(3));
Hfraction=(QFratio(2)+QFratio(3))/(QFratio(1)+QFratio(2)+QFratio(3));

%%%do simulation:
M=3000; %number of simulating molecules
N=END-START+1;
HDmatrix=ones(M,N); %0=H; 1=D
kmax=max([max(kcDH(START:END)), max(kcHD(START:END))]);
deltaT=0.1/kmax;   %step size of simulation time
for time=0:deltaT:Time_p
    ran1=rand(M,N);
    for i=1:M
        for j=START:END
            switch HDmatrix(i,j-START+1)
                case 0 %H
                    if ran1(i,j-START+1)<=(1-exp(-kcHD(j)*deltaT/pf(j)))*Dfraction
                        HDmatrix(i,j-START+1)=1; %H->D
                    end
                case 1 %D
                    if ran1(i,j-START+1)<=(1-exp(-kcDH(j)*deltaT/pf(j)))*Hfraction
                        HDmatrix(i,j-START+1)=0; %D->H
                    end
            end
        end
    end
    disp(time)
end

%%%consider N-term two residues and Prolines
for i=1:N
    if i<3 || Seq(i+START-1)=='P'
        HDmatrix(:,i)=zeros(M,1); %all must be H
    end
end

%%%get simulation result distribution
Distr=zeros(1,maxD+1);
deltaMass=zeros(1,M);
for i=1:M
    deltaMass(i)=sum(HDmatrix(i,:));
    Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Distr(x) means x-1 units of mass above monoisotopic
end
Distr=Distr/sum(Distr); %normalization

%%%do convolution with allH peaks:
obsDistr=conv(distND, Distr);
obsDistr=obsDistr/sum(obsDistr); %normalization

%%%get MS observable peaks:
clear obsPeaks
DM=1.00628; %delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
obsPeaks(1,1)=peptideMass/Charge+1.007276; %m/z of mono; 1.007276 is the mass of proton
obsPeaks(1,2)=obsDistr(1);
for i=2:(maxD+maxND+1)
    obsPeaks(i,1)=obsPeaks(i-1,1)+DM/Charge;
    obsPeaks(i,2)=obsDistr(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.axes1,'NextPlot','replace')
stem(handles.axes1,obsPeaks(:,1),obsPeaks(:,2),'k')


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
table1 = get(handles.uitable1,'Data');
Temp=table1{1,2};
pH=table1{3,2};

table2 = get(handles.uitable2,'Data');
clear Seq
for i=1:size(table2,1)
    Seq(i)=table2{i,2};
end
kcDH = fbmme_dh(Seq, pH, Temp, 'poly'); %call fbmme_dh.m
R=1.9858775; %gas constant(cal/(K*mol)); the value in spreadsheet is 1.987
for i=1:size(table2,1)
    table2{i,3} = kcDH(i); %update kcDH upon Temp or pH change
    table2{i,5} = R*(Temp+273.15)*log(table2{i,4})/1000; %update deltaG upon Temp change
end
set(handles.uitable2,'Data',table2);



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

table1 = get(handles.uitable1,'Data');
Temp=table1{1,2};

table2 = get(handles.uitable2,'Data');
R=1.9858775; %gas constant(cal/(K*mol)); the value in spreadsheet is 1.987
table2{eventdata.Indices(1),5} = R*(Temp+273.15)*log(table2{eventdata.Indices(1),4})/1000; %update deltaG
set(handles.uitable2,'Data',table2);



