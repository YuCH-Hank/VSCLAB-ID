clc, clear, close all;

%  2022.06.14, Hank Yu

global Show_Saved_File


%% Change Here ======================================================
Show_Saved_File = false;                                                     ... 是否顯示已儲存的檔案
% 讀取檔案之位置 -----------------------------------------------------------
Load_Path    = 'Text\';

% 儲存檔案之位置 -----------------------------------------------------------
Save_Path    = {'Text\Complex\'};

% 訂定變數名稱 (順序不可調整) -----------------------------------------------
% 系統阻尼, 系統剛性, 馬達阻尼, 負載阻尼, 馬達庫倫摩擦, 馬達轉動慣量, 齒數比 ---
SYM_FJ_Var   = {'C', 'K', 'Dl', 'Dm', 'fm', 'mM'}; 

% 訂定迴歸矩陣收集參數順序 (順序可以調整, 參數名稱需與自定義的變數相同) ---------
SYM_Collect_Order = {'C', 'K', 'mM', 'Dm', 'fm', 'Dl', 'Iner'};

% 迴歸矩陣 上下堆疊順序 -----------------------------------------------------
W_Order = {'m', 'l'};
    
% Flexible Joint Gear Ratio -----------------------------------------------
GearRatio = [1/50; 1/50];

% Flexible Joint State (Matlab 計算時所使用的系統狀態) ----------------------
SYM_FJ_State = {'pl', 'vl', 'al', 'pm', 'vm', 'am', 'NG'};

% Function Flexible Joint State (輸出成 Matlab function 所定義的系統變數)----
Function_FJ_State = {'Pl', 'Vl', 'Al', 'Pm', 'Vm', 'Am', 'NG', 'Beta'}; 

%% Step 1: Create Variables =========================================
% Expand Variables --------------------------------------------------------
ExpandStuff = {'ax', 'Beta', 'SYM_Geometric', 'SYM_State', 'SYM_SYMORO', 'Function_State', ...
               'W'};
load([Load_Path, 'rec.mat'])
for i = 1 : length(ExpandStuff)
    s = ExpandStuff{i};
    eval([s, '= rec.(s);'])
end

% Create State and SYMORO Variables ---------------------------------------
for i = 1 : ax
    for j = 1 : length(SYM_SYMORO)
        syms([SYM_SYMORO{j}, num2str(i)]);
    end
    for j = 1 : length(SYM_State)
        syms([SYM_State{j}, num2str(i)]);
    end
end

% Create Geometric Variables ----------------------------------------------
for i = 1 : size(SYM_Geometric, 1)
    eval(['syms ', SYM_Geometric{i, 1}])
end

% 紀錄變數順序 -------------------------------------------------------------
Collect_Order = zeros(1, length(SYM_Collect_Order));
for i = 1 : length(SYM_Collect_Order)
    z = 1 : length(SYM_Collect_Order);
    if strcmpi(SYM_Collect_Order{i}, 'Iner')
        Collect_Order(i) = nan;
    else
        Collect_Order(i) = z(strcmpi(SYM_FJ_Var, SYM_Collect_Order{i}));
    end
end

% 製作變數矩陣 -------------------------------------------------------------
SYM  = [SYM_FJ_Var(:)', SYM_FJ_State(:)'];
var_ = {'C', 'K', 'Dl', 'Dm', 'fm', 'mM', ...
        'pl', 'vl', 'al', 'pm', 'vm', 'am', 'NG'};

for j = 1 : length(SYM)
    tmp = [];
    for i = 1 : ax
        syms([SYM{j}, num2str(i)]);
        tmp = [tmp; eval([SYM{j}, num2str(i)])];
    end
    eval([SYM{j}, '= tmp;'])
    eval(['v_', var_{j}, '= tmp;'])
end

%% Step 2: Define Flexible Force ------------------------------------
tau_s = v_K .* (v_NG .* v_pm - v_pl) + v_C .* (v_NG .* v_vm - v_vl);
W_K_m = diag(v_NG .*(v_NG .* v_pm - v_pl ));
W_K_l = diag( - ( v_NG .* v_pm - v_pl ));

W_C_m = diag(v_NG .*(v_NG .* v_vm - v_vl ));
W_C_l = diag( - ( v_NG .* v_vm - v_vl ));

%% Step 3: Inverse Dynamic Model ====================================
%  把迴歸矩陣中有關於 FSi, FVi, IAi 項目 全部都刪掉 i = 1, ..., ax -----------

Index = true(ones(1, length(Beta)));
for i = 1 : ax
    Index = Index & ~strcmpi(Beta, [SYM_SYMORO{end}, num2str(i)]) ...
                  & ~strcmpi(Beta, [SYM_SYMORO{end - 1}, num2str(i)]) ...
                  & ~strcmpi(Beta, [SYM_SYMORO{end - 2}, num2str(i)]);
end

v_Iner = Beta(Index);
W_Iner_l = W(:, Index); clear W
W_Iner_m = zeros(size(W_Iner_l));

% 將慣量迴歸矩陣內的 p, v, a 置換成負載端 pl, vl, al -------------------------
for i = 1 : ax
    for j = 1 : 3
        W_Iner_l = subs( W_Iner_l, ...
                       str2sym([ SYM_State{j}, num2str(i) ]), ...
                       str2sym([ SYM_FJ_State{j}, num2str(i) ])) ;
    end
end

% 檢查 IAi 是否在變數內， 若無則取代 ZZRi 為 ZZRi - IAi ----------------------
for i = 1 : ax
    s_ax = num2str(i);
    s_IA = SYM_SYMORO{end - 2};
    if sum(strcmpi(Beta, [s_IA, s_ax])) == 0
        s_ZZR = SYM_SYMORO{12};
        v_Iner(strcmpi(v_Iner, [s_ZZR, s_ax])) = {[s_ZZR, s_ax, '-', s_IA, s_ax]};
    end
end

v_Iner = str2sym(v_Iner);

% 獲得慣量參數以及慣量迴歸矩陣 -----------------------------------------------
StrDisp('Inertia Parameters', '=')
disp(v_Iner)
StrDisp('W_Iner', ' ')
disp(W_Iner_l)

% 製作摩擦力, 馬達慣量之迴歸矩陣 ---------------------------------------------
W_Dl_m = zeros(ax);         W_Dl_l = diag(v_vl);
W_Dm_m = diag(v_vm);        W_Dm_l = zeros(ax);
W_fm_m = diag(sign(v_vm));  W_fm_l = zeros(ax);
W_mM_m = diag(v_am);        W_mM_l = zeros(ax);

% 重組迴歸矩陣 -------------------------------------------------------------
[W_Full_l, W_Full_m, Beta] = deal([]);
for i = 1 : length(Collect_Order)
    if isnan(Collect_Order(i))
        W_Full_l = [W_Full_l, eval('W_Iner_l')];
        W_Full_m = [W_Full_m, eval('W_Iner_m')];
        Beta   = [Beta,   eval('v_Iner')];
    else
        W_Full_l = [W_Full_l, eval(['W_', var_{Collect_Order(i)}, '_l'])];
        W_Full_m = [W_Full_m, eval(['W_', var_{Collect_Order(i)}, '_m'])];
        Beta   = [Beta,   eval(['v_', var_{Collect_Order(i)}, char(".'")])];
    end
end

W = eval(['[W_Full_', W_Order{1} '; W_Full_', W_Order{2}, '];']);

% 獲得慣量參數以及慣量迴歸矩陣 -----------------------------------------------
StrDisp('All Parameters', '=')
disp(Beta)
StrDisp('W_Full', ' ')
disp(W)

%% Step 4: Inverse Dynamic Model for Motor and Load Side ============
[W_l, W_m, Beta_l, Beta_m] = deal([]);
for i = 1 : length(Collect_Order)
    if isnan(Collect_Order(i))
        W_l    = [W_l, eval('W_Iner_l')];
        Beta_l = [Beta_l,   eval('v_Iner')];
    else
        if strcmpi(var_{Collect_Order(i)}, 'Dl')
            W_l    = [W_l,    eval(['W_', var_{Collect_Order(i)}, '_l'])];
            Beta_l = [Beta_l, eval(['v_', var_{Collect_Order(i)}, char(".'")])];
        else
            W_m    = [W_m,    eval(['W_', var_{Collect_Order(i)}, '_m'])];
            Beta_m = [Beta_m, eval(['v_', var_{Collect_Order(i)}, char(".'")])];
        end
    end
end

% 獲得慣量參數以及慣量迴歸矩陣 -----------------------------------------------
StrDisp('Motor Side Parameters', '=')
disp(Beta_l)
StrDisp('W_l', ' ')
disp(W_l)

StrDisp('Link Side Parameters', '=')
disp(Beta_m)
StrDisp('W_m', ' ')
disp(W_m)

%% Step 5: Validation With DDM ======================================
[M_l, N_l] = deal(rec.M, rec.N);

for i = 1 : ax
    for j = 1 : 3
        M_l = subs( M_l, str2sym([ SYM_State{j}, num2str(i) ]), ...
                         str2sym([ SYM_FJ_State{j}, num2str(i) ])) ;
        N_l = subs( N_l, str2sym([ SYM_State{j}, num2str(i) ]), ...
                         str2sym([ SYM_FJ_State{j}, num2str(i) ])) ;
    end
end

for i = 1 : ax
    v_IA(i, i) = eval([SYM_SYMORO{end - 2}, num2str(i)]);
end
M_l = M_l - v_IA;
M_m = diag(v_mM);
N_l = N_l + (v_Dl .* v_vl );
N_m = (v_Dm .* v_vm + v_fm .* sign(v_vm) + v_NG .* tau_s); 

M = eval(['[M_', W_Order{1} ', zeros(ax); zeros(ax), M_', W_Order{2}, '];']);
N = eval(['[N_', W_Order{1} '; N_', W_Order{2}, '];']);

% DDM v.s. IDM ------------------------------------------------------------
ddm = M * [v_am; v_al] + N - [zeros(ax, 1) ; tau_s];
idm = W * Beta .';
idm2 = [W_m * Beta_m .'; W_l * Beta_l .' - tau_s ;];

% Cheeck full Regression --------------------------------------------------
StrDisp(['Validate Full Regression'])
CheckDynamic(ddm, idm, ax);

StrDisp(['Validate Seperate Regression'])
CheckDynamic(ddm, idm2, ax);

%% Step 6. SaveFunctions ============================================
% 變數變換 P, V, A, Beta ---------------------------------------------------
stuff = {'W', 'W_m', 'W_l', 'M', 'N'};

FirstVariables = '~'; % '[]', 'this', ''                                   ... 是否需要建立第一項變數

% 資料蒐集 -----------------------------------------------------------------
CollectStuff = {'ax', 'Beta', 'SYM_Geometric', 'SYM_State', 'SYM_SYMORO', 'Function_State', 'FirstVariables'};
for i = 1 : length(CollectStuff)
    s = CollectStuff{i};
    rec.(s) = eval(s);
end

% Save Functions ----------------------------------------------------------
for i = 1 : length(Save_Path)
    disp(' ')
    % Create Folder -------------------------------------------------------
    fd = Save_Path{i};
    MakeDir(fd);
    
    % Create Beta.txt -----------------------------------------------------
    FileName = [fd, 'Beta_Flexible.txt'];
    Content  = Beta;
    CreateTxtFile(FileName, Content);
       
    % Create Each Functions -----------------------------------------------
    for j = 1 : length(stuff)
        s = stuff{j};
        InputName = [];
        switch s
            case {'W'} % Pl, Vl, Al, Pm, Vm, Am
                InputName = CreateInputName([1 : 6], Function_FJ_State);

            case {'W_m'} % Pl, Vl, Al, Pm, Vm, Am
                InputName = CreateInputName([1 : 2, 4:6], Function_FJ_State);
                
            case {'W_l'} % Pl, Vl, Al, Pm, Vm, Am
                InputName = CreateInputName([1 : 3], Function_FJ_State);
                
            case {'M'} % Pl, Pm, Beta
                InputName = CreateInputName([1, 8], Function_FJ_State);
                
            case 'N' % Pl, Vl, Pm, Vm, Beta
                InputName = CreateInputName([1, 2, 4, 5, 8], Function_FJ_State);
        end
        
        if ~isempty(FirstVariables)
            InputName = [FirstVariables, ' , ', InputName ];
        end
        
        rec.(s) = eval(s);
        
        % 變數變換 ---------------------------------------------------------
        for k = 1 : length(SYM_FJ_State)     % 不同的變換項目
            fun = @(s, i) subs(s, str2sym([SYM_FJ_State{k}, num2str(i)]), str2sym([Function_FJ_State{k}, '(', num2str(i), ')']) );
            for t = 1 : ax % 變換項目數量 ----------------------------------
                eval([s, '= fun(', s, ', t);']);
            end
        end
        
        fun = @(s, i) subs(s, Beta(i), str2sym([Function_FJ_State{end}, '(', num2str(i), ')']) );
        for t = 1 : length(Beta)
            eval([s, '= fun(', s, ', t);']);
        end
        
        % 資料儲存 ---------------------------------------------------------
        FunctionName = ['Flexible_', s, '_Full'];
        FileName     = [fd, FunctionName, '.m'];
        Content      = eval(s);
        CreateMFile(FileName, FunctionName, InputName, Content, SYM_Geometric, {Function_FJ_State{end - 1}, GearRatio});
        
    end
%     save([fd, 'rec.mat'], 'rec');
end

disp(' ')


%% Function =========================================================
function CreateTxtFile(FileName, Content)
global Show_Saved_File
File = fopen( FileName , 'wt' ) ;
PrintContent( File, Content )
fclose( File ) ;
if Show_Saved_File, open(FileName); end
StrDisp([FileName, ' has bean created'])
end

function CreateMFile(FileName, FunctionName, Input, Content, SYM_Geometric, GearRatio)
global Show_Saved_File
File = fopen( FileName , 'wt' ) ;
fprintf( File , ['function [ Output ] = ', FunctionName, '( ', Input, ' )\n']);

for i = 1 : size(SYM_Geometric, 1)
    fprintf( File , [SYM_Geometric{i, 1}, ' = ', num2str(SYM_Geometric{i, 2}), ';\n']) ;
end

fprintf( File , [GearRatio{1}, ' = [']) ;
tmp = num2str(GearRatio{2});
for i = 1 : size(tmp, 1)
    fprintf( File , ['[', tmp(i, :), ' ];\n']) ;
end
fprintf( File , '];\n') ;

fprintf( File , 'Output = [...\n');
PrintContent(File, Content);
fprintf( File , '];\n\n');
fprintf( File , 'end\n');
fclose( File ) ;

if Show_Saved_File, open(FileName); end
StrDisp([FileName, ' has bean created'])
end

function PrintContent(File, Content)
for i = 1 : size( Content , 1 )
    
    fprintf( File , '%s \t' , char( Content( i , : ) ) ) ;
    fprintf( File , '\n' ) ;
    
end
end

function StrDisp(str, varargin)
res = '>> -----------------------------------------------------------------------';
res(4 : 3 + length(str) + 1) = [str, ' '];
if ~isempty(varargin)
    if length(varargin) == 2 && ~varargin{2}
        res = [res(1 : 3 + length(str)), res(3 + length(str) + 2 : end), '-'];
    end
    res = strrep(res, '-', varargin{1});
end

disp(res);
end

function CheckDynamic(ddm, idm, ax)
[Bool_Motor, Bool_Link] = deal(false);
for Ax = 1 : ax
    Bool_Motor = Bool_Motor & (simplify(ddm(Ax, 1) - idm(Ax, 1)) == 0);
    Bool_Link  = Bool_Link  & (simplify(ddm(Ax + ax, 1) - idm(Ax + ax, 1)) == 0);
end

if ~Bool_Motor && ~Bool_Link
    StrDisp('OK')
else
    if Bool_Link
        StrDisp('Please Check Link Side');
    else
        StrDisp('Please Check Motor Side');
    end
end
end

function varargout = MakeDir(Path)
if ~exist(Path, 'dir')
    mkdir(Path)
    StrDisp([Path, ' Created Successfully'])
    varargout = {false};
else
    StrDisp([Path, ' Existed'])
    varargout = {true};
end
end

function InputName = CreateInputName(state, Function_State)
InputName = [];
for t = 1 : length(state)
    InputName = [InputName, Function_State{state(t)}];
    if t ~= length(state)
        InputName = [InputName, ' , '];
    end
end
end
