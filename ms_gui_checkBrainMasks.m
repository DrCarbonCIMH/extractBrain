function varargout = ms_gui_checkBrainMasks(varargin)
% MS_GUI_CHECKBRAINMASKS MATLAB code for ms_gui_checkBrainMasks.fig
%      MS_GUI_CHECKBRAINMASKS, by itself, creates a new MS_GUI_CHECKBRAINMASKS or raises the existing
%      singleton*.
%
%      H = MS_GUI_CHECKBRAINMASKS returns the handle to a new MS_GUI_CHECKBRAINMASKS or the handle to
%      the existing singleton*.
%
%      MS_GUI_CHECKBRAINMASKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MS_GUI_CHECKBRAINMASKS.M with the given input arguments.
%
%      MS_GUI_CHECKBRAINMASKS('Property','Value',...) creates a new MS_GUI_CHECKBRAINMASKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ms_gui_checkBrainMasks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ms_gui_checkBrainMasks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ms_gui_checkBrainMasks

% Last Modified by GUIDE v2.5 25-Jan-2018 06:02:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ms_gui_checkBrainMasks_OpeningFcn, ...
                   'gui_OutputFcn',  @ms_gui_checkBrainMasks_OutputFcn, ...
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


% --- Executes just before ms_gui_checkBrainMasks is made visible.
function ms_gui_checkBrainMasks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ms_gui_checkBrainMasks (see VARARGIN)

% Choose default command line output for ms_gui_checkBrainMasks
handles.output = hObject;
handles.spmCheckDone=0;
% we need spm12!!
if isempty(strfind(which('spm'),'spm12'));
    uiwait(warndlg('ms_gui_checkBrainMasks needs spm12! I will switch..', 'Switch to spm12'));
    spm12p; 
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ms_gui_checkBrainMasks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ms_gui_checkBrainMasks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function G_slider_Callback(hObject, eventdata, handles)
% hObject    handle to G_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.G_text.String=num2str(round(handles.G_slider.Value));
plotBrainMask(handles);

% --- Executes during object creation, after setting all properties.
function G_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Accept_pushbutton.
function Accept_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Accept_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mtx=handles.Mtx;
G = str2double(handles.G_text.String);
I_border = handles.I_border;

maskBrain=zeros(size(Mtx));
if iscell(I_border)
    for n=1:size(Mtx,3)
        maskBrain(:,:,n)=I_border{G}{n};
    end
else
    maskBrain = squeeze(I_border(G,:,:,:));
end
[d,name,ext]=fileparts(handles.BC);
Vtmp=spm_vol(handles.BC);
Vtmp.fname = [d filesep name '_brain' ext];
fprintf('writing %s\n', Vtmp.fname)
spm_write_vol(Vtmp, Mtx.*maskBrain);
Vtmp.fname = [d filesep name '_brainmask' ext];
fprintf('writing %s\n', Vtmp.fname)
spm_write_vol(Vtmp, maskBrain);
fprintf(' done\n');

% --- Executes on button press in Load_pushbutton.
function Load_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Load_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
P=spm_select(Inf,'mat','Select BrainMasks mat-files',[],pwd,'BrainMasks*');
handles.P=cellstr(P);
handles.selData_text.String='1';
handles=preparePlotting(handles);

handles.spmCheckDone=0;
plotBrainMask(handles);
handles.spmCheckDone=1;
guidata(hObject, handles);

function handles=preparePlotting(handles)
fprintf('Loading %s \n', handles.P{str2double(handles.selData_text.String)})
load(handles.P{str2double(handles.selData_text.String)});
handles.BC=BC;
[~,name,~]=fileparts(BC); handles.DataName_text.String=name;
handles.I_border=I_border;
handles.optG=optG;
delete(handles.axes1.Children)
idx=find(G_I>0,1,'last'); if idx~=length(G_I); endOfGI= idx-1; else; endOfGI=idx; end
G_I = G_I(1:endOfGI)/1000; handles.G_I = G_I; 
plot(G_I,'b', 'Parent', handles.axes1); plot(handles.axes1,optG,handles.G_I(optG),'rx', 'Tag','optG');
handles.G_slider.Min=1; handles.G_slider.Max=length(G_I); handles.G_slider.Value=optG;
handles.G_slider.SliderStep = [1/(handles.G_slider.Max-1) , 1/(handles.G_slider.Max-1) ];
handles.G_text.String=num2str(optG);
handles.Mtx=spm_read_vols(spm_vol([BC ',1']));
handles.OptG_text.String = sprintf('OptG: %.0f',optG);
V=spm_vol(BC); V=V(1); voxdim=spm_imatrix(V.mat); handles.voxdim=abs(voxdim(7:9));

function plotBrainMask(handles)
Mtx=handles.Mtx;
G = str2double(handles.G_text.String); % get the selected G
I_border = handles.I_border;

maskBrain=zeros(size(Mtx));
if iscell(I_border)
    for n=1:size(Mtx,3)
        maskBrain(:,:,n)=I_border{G}{n};
    end
else
    maskBrain = squeeze(I_border(G,:,:,:));
end
handles.vol_text.String = sprintf('Brain volume: %.2f mm^3\n', length(find(maskBrain))*prod(handles.voxdim)/1000);
% plot the GI
delete(findobj(handles.axes1,'Tag','G'))
plot(handles.axes1,G,handles.G_I(G),'ro', 'Tag','G');

% use spm_check_reg
[d,name,ext]=fileparts(handles.BC);
Vtmp=spm_vol(handles.BC); 
if length(Vtmp)>1; Vtmp=spm_vol([handles.BC ',1']); end
Vtmp.fname = [d filesep 'brainTMP' ext];
% spm_write_vol(Vtmp, Mtx.*maskBrain);
spm_write_vol(Vtmp, maskBrain);
% add the difference compared to the found G
maskBrainF=zeros(size(Mtx));
if iscell(I_border)
    for n=1:size(Mtx,3)
        maskBrainF(:,:,n)=I_border{handles.optG}{n};
    end
else
    maskBrainF = squeeze(I_border(handles.optG,:,:,:));
end
diffMask = abs(maskBrain-maskBrainF); diffMask(diffMask==0)=NaN;
Vtmp2=Vtmp; Vtmp2.fname = [d filesep 'brainTMP2' ext];
spm_write_vol(Vtmp2, Mtx.*diffMask);
% that is a pain in the a**!!!!
% toWrite(:,:,:,1)=maskBrain;
% toWrite(:,:,:,2)=diffMask;
% % Vtmp.dim=[Vtmp.dim size(toWrite,4)];
% for ix=1:size(toWrite,4)
%     Vtmp(ix)=Vtmp(1);
%     Vtmp(ix).private.dat.dim=size(toWrite);
%     Vtmp(ix).private.dat=toWrite;
%     spm_write_vol(Vtmp(ix), toWrite(:,:,:,ix));
% end

% show something
% global st; xhair=st.centre;
if ~handles.spmCheckDone || isempty(findobj(0,'type','figure','tag', 'Graphics'))
%     spm_check_registration(char([handles.BC ',1'], Vtmp.fname),{'Red: curr. sel. mask; yell: diff', 'yell: difference'});
    spm_check_registration(char([handles.BC ',1'], Vtmp.fname));
end
% spm_check_registration(char(handles.BC))
addColImage(Vtmp,Vtmp2, handles.trans_slider.Value);


function addColImage(V,V2, trans)
global st

st.vols{1}.mapping='histeq';

% V is the original BrainMask; V2 the difference
% st.vols{1}.blobs{1}.vol=V;
% st.vols{1}.blobs{1}.mat=V.mat;
% st.vols{1}.blobs{1}.max=1;
% st.vols{1}.blobs{1}.min=0.01;
% st.vols{1}.blobs{1}.colour = struct('cmap',[0 0 0; 1 0 0],'prop',trans);

% st.vols{1}.blobs{2}.vol=V2;
% st.vols{1}.blobs{2}.mat=V2.mat;
% st.vols{1}.blobs{2}.max=1;
% st.vols{1}.blobs{2}.min=0;
% st.vols{1}.blobs{2}.colour = struct('cmap',[0 0 0; 1 1 0],'prop',trans);

% st.vols{2}.blobs{1}=st.vols{1}.blobs{2};
% % st.vols{2}.blobs{1}.colour = struct('cmap',[0 0 0; 1 1 0],'prop',0.3);

% st.vols{2}.blobs{1}.vol=V2;
% st.vols{2}.blobs{1}.mat=V2.mat;
% st.vols{2}.blobs{1}.max=1;
% st.vols{2}.blobs{1}.min=0;
% st.vols{2}.blobs{1}.colour = [1 1 0];

% for contour
st.vols{2}.contour.images=1; % to initialise that there should be a contour of image 2 in image 1
% for d=1:3
% CData=st.vols{2}.ax{d}.d.CData;
% CData(CData>0)=1; st.vols{2}.ax{d}.d.CData=CData;
% end

hM = findobj(st.vols{2}.ax{1}.cm,'Label','Contour');
hM.UserData.nblines=1;
spm_ov_contour('redraw',2,1)
% delete(hM)
% contour_display(2,NaN)
% force a "reposition"
spm_orthviews('reposition');
% contour_delete(1);
% st.vols{1}.blobs{1}.cbar

function contour_display(i,o)

global st

if nargin < 2
    o = spm_input('Select image(s)', '!+1', 'e', ...
        num2str(spm_orthviews('valid_handles')));
    o = intersect(spm_orthviews('valid_handles'),o);
elseif isinf(o)
    o = spm_orthviews('valid_handles');
elseif isnan(o)
    o = setxor(spm_orthviews('valid_handles'),i);
end

% try
%     hM = findobj(st.vols{i}.ax{1}.cm,'Label','Contour');
%     UD = get(hM,'UserData');
%     nblines = UD.nblines;
%     linestyle = UD.style;
% catch
%     nblines = 3;
%     linestyle = 'r-';
% end
linewidth = 1;
nblines = 1;
linestyle = 'r-';
contour_delete(i);

lh = {};
sw = warning('off','MATLAB:contour:ConstantData');
for d = 1:3
    CData = sqrt(sum(get(st.vols{i}.ax{d}.d,'CData').^2, 3));
    CData(isinf(CData)) = NaN;
    CData(isnan(CData)) = 0;
    CData(CData>0) = 1;
    for h = o(:)'
        set(st.vols{h}.ax{d}.ax,'NextPlot','add');
        [C,lh{end+1}] = ...
            contour(st.vols{h}.ax{d}.ax,CData,...
            nblines,linestyle,'LineWidth',linewidth);
    end
end
warning(sw);
set(cat(1,lh{:}),'HitTest','off');

st.vols{i}.contour.images = o;
st.vols{i}.contour.handles = lh;

%==========================================================================
function contour_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD

global st

contour_delete(i);
contour_display(i,st.vols{i}.contour.images);


%==========================================================================
function contour_delete(i)

global st

try, delete(cat(1,st.vols{i}.contour.handles{:})); end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'BC') 
    [d,name,ext]=fileparts(handles.BC);
    delete([d filesep 'brainTMP' ext]); 
    delete([d filesep 'brainTMP2' ext]);
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in plus_pushbutton.
function plus_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plus_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'BC') 
    [d,name,ext]=fileparts(handles.BC);
    delete([d filesep 'brainTMP' ext]); 
    delete([d filesep 'brainTMP2' ext]);
end
newNum=str2double(handles.selData_text.String)+1;
handles.selData_text.String=num2str(newNum);
handles=preparePlotting(handles);
handles.spmCheckDone=0;
plotBrainMask(handles);
handles.spmCheckDone=1;
guidata(hObject, handles);

% --- Executes on button press in minus_pushbutton.
function minus_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to minus_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'BC') 
    [d,name,ext]=fileparts(handles.BC);
    delete([d filesep 'brainTMP' ext]); 
    delete([d filesep 'brainTMP2' ext]);
end
newNum=str2double(handles.selData_text.String)-1;
handles.selData_text.String=num2str(newNum);
handles=preparePlotting(handles);
handles.spmCheckDone=0;
plotBrainMask(handles);
handles.spmCheckDone=1;
guidata(hObject, handles);


% --- Executes on slider movement.
function trans_slider_Callback(hObject, eventdata, handles)
% hObject    handle to trans_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global st
trans = get(hObject,'Value');
st.vols{1}.blobs{1}.colour = struct('cmap',[0 0 0; 1 0 0],'prop',trans);
st.vols{1}.blobs{2}.colour = struct('cmap',[0 0 0; 1 1 0],'prop',trans);
spm_orthviews('reposition');

% --- Executes during object creation, after setting all properties.
function trans_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
