function varargout = digital_filters(varargin)
% DIGITAL_FILTERS MATLAB code for digital_filters.fig
%      DIGITAL_FILTERS, by itself, creates a new DIGITAL_FILTERS or raises the existing
%      singleton*.
%
%      H = DIGITAL_FILTERS returns the handle to a new DIGITAL_FILTERS or the handle to
%      the existing singleton*.
%
%      DIGITAL_FILTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGITAL_FILTERS.M with the given input arguments.
%
%      DIGITAL_FILTERS('Property','Value',...) creates a new DIGITAL_FILTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before digital_filters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to digital_filters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help digital_filters


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @digital_filters_OpeningFcn, ...
                   'gui_OutputFcn',  @digital_filters_OutputFcn, ...
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


% --- Executes just before digital_filters is made visible.
function digital_filters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to digital_filters (see VARARGIN)


%draw unit circle
c_draw(hObject, eventdata, handles);

 %Saw tooth
% T = 10*(1/50);
% st_fs = 1000;
% dt = 1/st_fs;
% st_t = 0:dt:T-dt;
% y = sawtooth(2*pi*50*st_t);
% handles.y_fft=fft(st,256);

%ECG signal
% [tm, y]=rdsamp('E:\Fantasia Database\f1o01.hea',2,3749,2500,1,0); %import data
global signal 
% sine wave
t=[0:0.1:20];
A=0.5;
f=1;
global flag
if flag ~= 1
    y=A*cos(f*t);
else 
    y = signal;
end
%calculate signal fft
handles.y_fft=fft(y);

%plot the original signal in time domain
axes(handles.axes4)
plot(y)
xlabel('Signal in time domain before filteration')


handles.p=[];   % Array holds poles in complex formula
handles.z=[];   % Array holds zeros in complex formula

handles.A=[];   % Array holds poles points in the listBox in string formula
handles.B=[];   % Array holds zeros points in the listBox in string formula
global choose
if choose ~= 2
    choose = 1
end
% Choose default command line output for digital_filters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes digital_filters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = digital_filters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PP.
function PP_Callback(hObject, eventdata, handles)
% hObject    handle to PP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global choose
global x y
%click to add a pole point in the unit circle
[x,y]=ginput(1);

%push a point and its conjugate to the poles array 
handles.p(length(handles.p)+1)=x+1j*y;
if choose == 2 
    handles.p(length(handles.p)+1)=x+1j*(-y);
end 
% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%push a point and its conjugate to the poles listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
if choose == 2
    handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];
end

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ZZ.
function ZZ_Callback(hObject, eventdata, handles)
% hObject    handle to ZZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global choose
global xzero yzero
%click to add a zero point in the unit circle
[xzero,yzero]=ginput(1);

%push a point and its conjugate to the zeros array 
handles.z(length(handles.z)+1)=xzero+1j*yzero;
if choose == 2
    handles.z(length(handles.z)+1)=xzero+1j*(-yzero);
end
% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%push a point and its conjugate to the zeros listbox array 
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
if choose == 2    
    handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(-yzero),')'];
end 
%show the point in the last point in the edit line

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in nP.
function nP_Callback(hObject, eventdata, handles)
% hObject    handle to nP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to remove a pole point in the unit circle
[x,y]=ginput(1);

%search for the selected pole from the poles array and remove it
temp=find(real(handles.p <(x+0.1)) & real(handles.p>(x-0.1)) );
handles.p(find(real(handles.p <(x+0.1)) & real(handles.p>(x-0.1)) ))=[];

%remove the selected pole from  the listBox
handles.A(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox1,'String',handles.A);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in nZ.
function nZ_Callback(hObject, eventdata, handles)
% hObject    handle to nZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to remove a zero point in the unit circle
[x,y]=ginput(1);

%search for the selected zero from the zero array and remove it
temp=find((real(handles.z <(x+0.1)) & real(handles.z>(x-0.1))));
handles.z(find((real(handles.z <(x+0.1)) & real(handles.z>(x-0.1)))))=[];

%remove the selected zero from  the listBox
handles.B(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing points
set(handles.listbox2,'String',handles.B);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear all arrays
handles.p=[];
handles.z=[];
handles.A=[];
handles.B=[];

%clear all plots
cla(handles.axes1,'reset')
cla(handles.axes2,'reset')
cla(handles.axes5,'reset')
cla(handles.axes9,'reset')

c_draw(hObject, eventdata, handles)

%update the listBoxes after clearing all
set(handles.listbox1,'String',handles.p);
set(handles.listbox2,'String',handles.z);
% Update handles structure
guidata(hObject, handles);



function c_draw(hObject, eventdata, handles)
%draw unit circle
circle_1 = exp(1i*(0:63)*2*pi/64); 
 axes(handles.axes1)
plot(real(circle_1),imag(circle_1),'.');

axis([-2 2 -2 2]); 
axis('equal'); 
hold on
plot( [0 0], [1.5 -1.5], '-')
plot( [1.5 -1.5], [0 0], '-')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
hold off;
% Update handles structure
guidata(hObject, handles);

function freq_plot(hObject, eventdata, handles)

%clear the unit circuit axes
cla(handles.axes1,'reset');
axes(handles.axes1)
c_draw(hObject, eventdata, handles);
hold on

%plot poles and zeros markers
plot_p=plot(real(handles.p),imag(handles.p),'X');
plot_z=plot(real(handles.z),imag(handles.z),'O');
set(plot_p,'markersize',8,'linewidth',2);
set(plot_z,'markersize',8,'linewidth',2);
hold off;

%Get the transfer function coeffecients
global a 
global b
[b,a]=zp2tf(handles.z',handles.p,1);

%Get the frequency response 
[h,w] = freqz(b,a,length(handles.y_fft));

%plot the frequency response mag 
axes(handles.axes2)
plot(w/pi,20*log10(abs(h)));
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
grid on;

%apply the filter in the orignial signal
%filter=h'.*handles.y_fft;
global signal
filt = filter(b,a,signal);
axes(handles.axes5)
plot(real(ifft(filt)))
xlabel('Signal in time domain after filteration')


unitcircle = exp(1i*(0:63)*pi/64)


%hena dah el manual gain
distance=ones(315,1);
for ii=1:length(unitcircle)
    for j=1:length(handles.z)        
      distance(ii)= distance(ii)*norm(unitcircle(ii)- handles.z(j));  
    end
end
 
for ii=1:length(unitcircle)
    for j=1:length(handles.p)       
      distance(ii)= distance(ii)*(1./norm(unitcircle(ii)- handles.p(j)));
    end
end
manual_gain=20*log(abs(distance));
axes(handles.axes9)  
plot(manual_gain);
xlim([0 60])
xticklabels({'0','0.2','0.4','0.5','0.6','0.8','1'})
yticklabels({''})



  
% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global choose 
  switch get(handles.popupmenu2,'Value')   
    case 1
     choose = 1;
    case 2
      choose = 2;
    otherwise
 end 
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global signal
global flag 
flag = 1
[filename,pathname]= uigetfile ( {'*.wav'} , 'Select File' ) ; 
 ExPath = fullfile(pathname, filename);
 [signal,fs] = audioread([pathname filename]);
 signal = signal(:,1);
 axes(handles.axes4)
 plot(signal)


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.


% --- Executes on button press in online.
function online_Callback(hObject, eventdata, handles)
% hObject    handle to online (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global signal
global a
global b
y = []
signal = signal(4000:end) ;
for n = 1:length(signal)
    Y = 0
    X = 0
    output = 0
    for j = 1:length(b)
        if n-j < 1
            break;
        end
        X = X + b(j) * signal(n-j); 
    end
    for k = 1:length(a)
        if n - k < 1
            break;
        end
        Y = a(k) * y(n-k);
    end
    output = (X - Y);
    y(end +1) = output;
    axes (handles.axes12)
    plot(y)
    drawnow
end

% Hint: get(hObject,'Value') returns toggle state of online


% --- Executes on button press in upp.
function upp_Callback(hObject, eventdata, handles)
global x y 
%search for the selected pole from the poles array and remove it
index = get(handles.listbox1 , 'value');
ss = strsplit(string(handles.A(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
x = str2double(ssx(2:end))
y = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.p)*100) == floor(str2num(ss)*100));

handles.p(temp)=[];
%remove the selected pole from  the listBox
handles.A(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox1,'String',handles.A);
if y > 1 || y <-1 
    warndlg('unstable filter ! ' , 'warning messsage ' ) 
    pause (1) ;
end
y=y+0.1 ; 
global choose 
if choose == 2
    %push a point and its conjugate to the poles array 
handles.p(length(handles.p)+1)=x+1j*y;
handles.p(length(handles.p)+1)=x+1j*(-y);
elseif choose == 1
    %push a point and its conjugate to the poles array 
    handles.p(length(handles.p)+1)=x+1j*y;
end


% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];
elseif choose == 1
    handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
    end

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);
% --- Executes on button press in down.

function downp_Callback(hObject, eventdata, handles)
global x y 

%search for the selected pole from the poles array and remove it
index = get(handles.listbox1 , 'value');
ss = strsplit(string(handles.A(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
x = str2double(ssx(2:end))
y = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.p)*100) == floor(str2num(ss)*100));

handles.p(temp)=[];
%remove the selected pole from  the listBox
handles.A(temp)=[];
% Update handles structure
guidata(hObject, handles);
%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
%show the listbox after removing  points
set(handles.listbox1,'String',handles.A);
if y > 1 || y <-1 
    warndlg('unstable filter ! ' , 'warning messsage ' ) 
    pause (1) ;
end
y=y-0.1 ; 
global choose 
if choose == 2 
    %push a point and its conjugate to the poles array 
handles.p(length(handles.p)+1)=x+1j*y;
handles.p(length(handles.p)+1)=x+1j*(-y);
elseif choose == 1
    %push a point and its conjugate to the poles array 
    handles.p(length(handles.p)+1)=x+1j*y;
end

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];
elseif choose == 1 
    handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
end

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);
% --- Executes on button press in rightp.
function rightp_Callback(hObject, eventdata, handles)
global x y 

%search for the selected pole from the poles array and remove it
index = get(handles.listbox1 , 'value');
ss = strsplit(string(handles.A(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
x = str2double(ssx(2:end))
y = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.p)*100) == floor(str2num(ss)*100));

handles.p(temp)=[];

%remove the selected pole from  the listBox
handles.A(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox1,'String',handles.A);
if x > 1 || x <-1 
    warndlg('unstable filter ! ' , 'warning messsage ' ) 
    pause (1) ;
end
x=x+0.1 ; 
global choose 
if choose == 2
    %push a point and its conjugate to the poles array 
handles.p(length(handles.p)+1)=x+1j*y;
handles.p(length(handles.p)+1)=x+1j*(-y);
elseif choose == 1
    %push a point and its conjugate to the poles array 
    handles.p(length(handles.p)+1)=x+1j*y;
end

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];
elseif choose == 1
    handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
end

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in leftp.
function leftp_Callback(hObject, eventdata, handles)
global x y 

%search for the selected pole from the poles array and remove it
index = get(handles.listbox1 , 'value');
ss = strsplit(string(handles.A(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
x = str2double(ssx(2:end))
y = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.p)*100) == floor(str2num(ss)*100));

handles.p(temp)=[];

%remove the selected pole from  the listBox
handles.A(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox1,'String',handles.A);
if x > 1 || x <-1 
    warndlg('unstable filter ! ' , 'warning messsage ' ) 
    pause (1) ;
end
x=x-0.1 ; 
global choose 
if choose == 2
    %push a point and its conjugate to the poles array 
handles.p(length(handles.p)+1)=x+1j*y;
handles.p(length(handles.p)+1)=x+1j*(-y);
elseif choose == 1
    %push a point and its conjugate to the poles array 
    handles.p(length(handles.p)+1)=x+1j*y;
end


% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];
elseif choose == 1
    handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
end

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in upz.
function upz_Callback(hObject, eventdata, handles)
% hObject    handle to upz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xzero yzero

%search for the selected pole from the poles array and remove it
index = get(handles.listbox2 , 'value');
ss = strsplit(string(handles.B(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
xzero = str2double(ssx(2:end))
yzero = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.z)*100) == floor(str2num(ss)*100));

handles.z(temp)=[];

%remove the selected pole from  the listBox
handles.B(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox2,'String',handles.B);

yzero=yzero+0.1 ; 
global choose 
if choose == 2
    %push a point and its conjugate to the zeros array 
handles.z(length(handles.z)+1)=xzero+1j*yzero;
handles.z(length(handles.z)+1)=xzero+1j*(-yzero);
elseif choose == 1
    %push a point and its conjugate to the zeros array 
    handles.z(length(handles.z)+1)=xzero+1j*yzero;
end


% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2 
%push a point and its conjugate to the poles listbox array 
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(-yzero),')'];
elseif choose == 1 
    handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
end

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in downz.
function downz_Callback(hObject, eventdata, handles)
global xzero yzero
%search for the selected pole from the poles array and remove it
index = get(handles.listbox2 , 'value');
ss = strsplit(string(handles.B(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
xzero = str2double(ssx(2:end))
yzero = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.z)*100) == floor(str2num(ss)*100));


handles.z(temp)=[];
%remove the selected pole from  the listBox
handles.B(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox2,'String',handles.B);

yzero=yzero-0.1 ; 
global choose 
if choose == 2
    %push a point and its conjugate to the zeros array 
handles.z(length(handles.z)+1)=xzero+1j*yzero;
handles.z(length(handles.z)+1)=xzero+1j*(-yzero);
elseif choose == 1
    %push a point and its conjugate to the zeros array 
    handles.z(length(handles.z)+1)=xzero+1j*yzero;
end


% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(-yzero),')'];
elseif choose == 1
    handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
end

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rightz.
function rightz_Callback(hObject, eventdata, handles)
global xzero yzero 

%search for the selected pole from the poles array and remove it
index = get(handles.listbox2 , 'value');
ss = strsplit(string(handles.B(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
xzero = str2double(ssx(2:end))
yzero = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.z)*100) == floor(str2num(ss)*100));


handles.z(temp)=[];

%remove the selected pole from  the listBox
handles.B(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox2,'String',handles.B);

xzero=xzero+0.1 ; 
global choose 
if choose == 2 
    %push a point and its conjugate to the zeros array 
handles.z(length(handles.z)+1)=xzero+1j*yzero;
handles.z(length(handles.z)+1)=xzero+1j*(-yzero);
elseif choose == 1
    %push a point and its conjugate to the zeros array 
    handles.z(length(handles.z)+1)=xzero+1j*yzero;
end


% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(-yzero),')'];
elseif choose == 1
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
end

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in leftz.
function leftz_Callback(hObject, eventdata, handles)
global xzero yzero 

%search for the selected pole from the poles array and remove it
index = get(handles.listbox2 , 'value');
ss = strsplit(string(handles.B(index)) , ',');
ssx = char(ss(1));
ssy = char(ss(2));
xzero = str2double(ssx(2:end))
yzero = str2double(ssy(1:end-1))
ss = char(ss(1));
ss = ss(2:end);
temp=find(floor(real(handles.z)*100) == floor(str2num(ss)*100));


handles.z(temp)=[];

%remove the selected pole from  the listBox
handles.B(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox2,'String',handles.B);

xzero=xzero-0.1 ; 
global choose 
if choose == 2 
    %push a point and its conjugate to the zeros array 
handles.z(length(handles.z)+1)=xzero+1j*yzero;
handles.z(length(handles.z)+1)=xzero+1j*(-yzero);
elseif choose == 1
    %push a point and its conjugate to the zeros array 
    handles.z(length(handles.z)+1)=xzero+1j*yzero;
end


% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
if choose == 2
%push a point and its conjugate to the poles listbox array 
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(-yzero),')'];
elseif choose == 1
handles.B{length(handles.B)+1}=['(',num2str(xzero),' , ',num2str(yzero),')'];
end

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);
