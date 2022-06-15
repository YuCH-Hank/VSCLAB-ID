clc; clear; close all;
%% SYMORO -> Matlab function ========================================
%  2022.06.14, Hank Yu

global Show_Saved_File

%% Change Here ======================================================
%  需要三個 SYMORO 檔案 XXX_dim.txt, XXX_base_inm.txt, XXX_regp.txt --------
File_Path       = 'Data\';                                                  ... SYMORO 檔案儲存位置 
Save_Path       = {'Text\', 'Text\Complex\'};                               ... 儲存檔案之位置
Show_Saved_File = false;                                                     ... 是否顯示已儲存的檔案
ax              = 2;                                                        ... 關節數量

% Select Different Functions 1 / 0 ----------------------------------------
Options = {'W', 'M', 'Mdiag', 'Moffdiag', 'Mdot', 'C', 'N', 'G', 'F';            ... 各函式是否需要計算
             1,   1,       1,          1,      1,   1,   1,   1,  1;};
         
% Geometric Variables -----------------------------------------------------
SYM_Geometric = {'D1', 0.24;
                 'D2', 0.24;
                 'g', -9.81};

% SYMORO Variables --------------------------------------------------------
SYM_State   = {'p', 'v', 'a'};                                              ... 位置 速度 加速度
SYM_SYMORO  = {'LXX', 'LXY', 'LXZ', 'LYY', 'LYZ', 'LZZ',                    ... 慣量參數
               'XXR', 'XYR', 'XZR', 'YYR', 'YZR', 'ZZR',                    ... 耦合慣量參數
               'MXR', 'MX', 'MYR', 'MY', 'MZR', 'MZ', 'M',                  ... 質量參數
               'IA', 'FV', 'FC'};                                           ... 馬達慣量、摩擦力
    
% Collect States ----------------------------------------------------------
Function_State = {'P', 'V', 'A', 'Beta'};                                   ... 儲存檔案的狀態變數
FirstVariables = '~'; % '[]', 'this', ''                                   ... 是否需要建立第一項變數

%% Step 1: Create Variables =========================================
% Create State and SYMORO Variables ---------------------------------------
for i = 1 : ax
    for j = 1 : length(SYM_SYMORO)
        syms([SYM_SYMORO{j}, num2str(i)]);
    end
    for j = 1 : length(SYM_State)
        syms([SYM_State{j}, num2str(i)]);
    end
end

% Create Options Switchs --------------------------------------------------
for i = 1 : length(Options)
    options.(Options{1, i}) = logical(Options{2, i});
end

% Create Geometric Variables ----------------------------------------------
for i = 1 : size(SYM_Geometric, 1)
    eval(['syms ', SYM_Geometric{i, 1}])
end

%% Step 2: IDM_Standard, 逆向動力學 ( XXX_dim.txt ) ==================
% Parse Data in files -----------------------------------------------------
File_data = dir([File_Path, '*_dim.txt']);
path = [File_data.folder, '\', File_data.name];

StrDisp('Start collecting Inverse Dynamic Model', '=')
StrDisp(['Load file from: ', File_data.name])
String = Parse(path);
try eval(String), catch, open(path); error(['Something Wrong in ', File_data.name]), end

% Check variables exist or not --------------------------------------------
Beta = {}; index = 1;
for i = 1 : ax
    for j = 1 : length(SYM_SYMORO) - 2
        for k = 1 : ax % DG
            name = ['DG', num2str(k), SYM_SYMORO{j}, num2str(i)];
            if exist(name, 'var') == 1
                Beta{index} = [SYM_SYMORO{j}, num2str(i)];
                index = index + 1;
                break;
            end
        end
    end
end

% Create Friction Variables behind ----------------------------------------
for j = 1 : 2
    for i = 1 : ax
        Beta{index} = [SYM_SYMORO{end - j + 1}, num2str(i)];
        index = index + 1;
    end
end

% show exist variables ----------------------------------------------------
StrDisp('Existed Variables')
disp(Beta)
StrDisp([num2str(length(Beta)), ' variables'], ' ')

if options.W % Create IDM -------------------------------------------------    
    StrDisp('Create IDM as: ', ' ')
    W = [];
    for i = 1 : ax % DG
        V = [];
        fprintf('   ');
        for j = 1 : length(Beta)
            name = ['DG', num2str(i), Beta{j}];
            if exist(name, 'var')
                V = [V, eval(name)];
                fprintf('%s', changelength(name, 9));
            else
                V = [V, 0];
                fprintf('%s', changelength('0', 9));
            end
        end
        fprintf('\n');
        W = [W; V];
    end
end

disp(' ')

%% Step 3. DDM_Base, 順向動力學 ( XXX_base_inm.txt ) =================
% Parse Data in files -----------------------------------------------------
File_data = dir([File_Path, '*_base_inm.txt']);
path = [File_data.folder, '\', File_data.name];

String = Parse(path);
try eval(String), catch, open(path); error(['Something Wrong in ', File_data.name]), end

% M -----------------------------------------------------------------------
M = [];
for i = 1 : ax
    rows = [];
    for j = 1 : ax 
        if j < i
            ind = [i, j];
        else
            ind = [j, i];
        end
        rows = [rows, eval(['A', num2str(ind(1)), num2str(ind(2))])];
    end
    M = [M; rows];
end

if options.Mdiag 
    Mdiag = diag( M ) ; 
end
if options.Moffdiag 
	Moffdiag = M - diag( diag( M ) ) ;  
end
if options.Mdot 
    Mdot = sym( zeros( ax , ax ) ) ;
    for i = 1 : ax
        for j = 1 : ax
            for k = 1 : ax
                Mdot( i , j ) = Mdot( i , j ) + diff( M( i , j ) , [ 'p' , num2str( k ) ] ) * [ 'v' , num2str( k ) ] ;
            end    
        end
    end
    Mdot = simplify( Mdot ) ;
end
if options.C 
    C = sym( zeros( ax , ax ) ) ;
    for i = 1 : ax
        for j = 1 : ax
            for k = 1 : ax

                C( i , j ) = C( i , j ) + ...
                    0.5 * (   diff( M( i , j ) , [ 'p' , num2str( k ) ] ) ... 
                            + diff( M( i , k ) , [ 'p' , num2str( j ) ] ) ...
                            - diff( M( j , k ) , [ 'p' , num2str( i ) ] ) ) * [ 'v' , num2str( k ) ] ;

            end    
        end
    end
    C = simplify( C ) ;
end
if options.N || options.F || option.G 
    File_data = dir([File_Path, '*_base_ccg.txt']);
    path = [File_data.folder, '\', File_data.name];
    String = Parse(path);
    try eval(String), catch, open(path); error(['Something Wrong in ', File_data.name]), end

    if options.N || option.G
        N = [];
        for i = 1 : ax
            N = [N; eval(['N', num2str(3), num2str(i)])];
        end
    end
    if options.G
        G = N ;
        for i = 1 : ax
            G = subs( G , [ SYM_State{2} , num2str( i ) ] , 0 ) ;
        end
    end
    if options.F
        F = [];
        for i = 1 : ax
            F = [F; eval([SYM_SYMORO{end}, num2str(i)]) * sign(eval(['v', num2str(i)])) + eval([SYM_SYMORO{end - 1}, num2str(i)]) * eval(['v', num2str(i)])];
        end
    end
end

%% Step 4: 收集基本動態參數 ( XXX_regp.txt ) ==========================
% Parse Data in files -----------------------------------------------------
File_data = dir([File_Path, '*_regp.txt']);
path = [File_data.folder, '\', File_data.name];

StrDisp('Start collecting Dynamic Parameters', '=')
StrDisp(['Load file from: ', File_data.name])
String = Parse(path);
try eval(String), catch, open(path); error(['Something Wrong in ', File_data.name]), end

% Collect Beta Full -------------------------------------------------------
BetaRegp  = cell2sym(Beta).';
Beta_Full = [];
for i = 1 : length(Beta)
   Beta_Full = [Beta_Full; simplify(eval(Beta{i}))]; 
end

StrDisp('Dynamic Parameter defined as below: ')
DataDisp({string(BetaRegp), string(Beta_Full)}, {'Beta', 'Beta_Full'})
disp(' ')

%% Step 5. Save Functions ===========================================
% 變數變換 P, V, A, Beta ---------------------------------------------------
changestuff = {ax,               @(s, i) subs(s, str2sym([SYM_State{1}, num2str(i)]), str2sym([Function_State{1}, '(', num2str(i), ')']) );
               ax,               @(s, i) subs(s, str2sym([SYM_State{2}, num2str(i)]), str2sym([Function_State{2}, '(', num2str(i), ')']) );
               ax,               @(s, i) subs(s, str2sym([SYM_State{3}, num2str(i)]), str2sym([Function_State{3}, '(', num2str(i), ')']) );
               length(BetaRegp), @(s, i) subs(s, BetaRegp(i),                         str2sym([Function_State{4}, '(', num2str(i), ')']) )};

stuff = {Options{1, logical([Options{2, :}])}};

% Save Functions ----------------------------------------------------------
for i = 1 : length(Save_Path)
    disp(' ')
    % Create Folder -------------------------------------------------------
    fd = Save_Path{i};
    MakeDir(fd);
    
    % Create Beta.txt -----------------------------------------------------
    FileName = [fd, 'Beta.txt'];
    Content  = BetaRegp;
    CreateTxtFile(FileName, Content);
    
    % Create Beta_Full.txt ------------------------------------------------
    FileName = [fd, 'Beta_Full.txt'];
    Content  = [BetaRegp, Beta_Full];
    CreateTxtFile(FileName, Content);
       
    % Create Each Functions -----------------------------------------------
    for j = 1 : length(stuff)
        s = stuff{j};
        
        switch s
            case {'M', 'Mdiag', 'Moffdiag', 'Mdot', 'G'} % P, Beta
                InputName = [Function_State{1}, ' , ', Function_State{4}];
                
            case {'C', 'N'} % P, V, Beta
                InputName = [Function_State{1}, ' , ', Function_State{2}, ' , ', Function_State{4}];
                
            case 'F' % V, Beta
                InputName = [Function_State{2}, ' , ', Function_State{4}];
                
            case 'W' % P, V, A
                InputName = [Function_State{1}, ' , ', Function_State{2}, ' , ', Function_State{3}];
        end
        
        if ~isempty(FirstVariables)
            InputName = [FirstVariables, ' , ', InputName ];
        end
        
        % 變數變換 ---------------------------------------------------------
        for k = 1 : length(changestuff)     % 不同的變換項目
            fun = changestuff{k, 2};

            for t = 1 : changestuff{k, 1}   % 變換項目數量
                eval([s, '= fun(', s, ', t);']);
            end
        end
        
        % 資料儲存 ---------------------------------------------------------
        FileName     = [fd, s, '_Full.m'];
        Content      = eval(s);
        FunctionName = [s, '_Full'];
        CreateMFile(FileName, FunctionName, InputName, Content, SYM_Geometric)
    end
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

function CreateMFile(FileName, FunctionName, Input, Content, SYM_Geometric)
global Show_Saved_File
File = fopen( FileName , 'wt' ) ;
fprintf( File , ['function [ Output ] = ', FunctionName, '( ', Input, ' )\n']);

for i = 1 : size(SYM_Geometric, 1)
    fprintf( File , [SYM_Geometric{i, 1}, ' = ', num2str(SYM_Geometric{i, 2}), ';\n']) ;
end

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

function new_name = changelength(name, LEN)
tmp = '                                               ';
new_name = tmp(1 : max(LEN, length(name)));
new_name(1:length(name)) = name;
end

function DataDisp(data, name)
tmp = data{2};
cols = length(name); rows = size(tmp, 1);
DISP = cell(rows, cols);

for i = 1 : cols
    tmp = reshape(data{i}, [rows, 1]);
    for j = 1 : rows
        DISP{j, i} = tmp(j, 1);
    end
end

table = cell2table(DISP);

table.Properties.VariableNames = [name(:)'];
disp(table);
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

function String = Parse(path)
fid = fopen(path);
tline = fgetl(fid);

String = []; start = false;
while ischar(tline)
    try
        tline = fgetl(fid); % Get Line
        if strcmpi(tline, 'Equations:')
            start = true;
        elseif strcmpi(tline, '*=*')
            start = false;
            break
        end
        if start && ~strcmpi(tline, 'Equations:')
            String = [String, tline];
        end
        
    catch
        disp(path)
        disp(tline)
        break;
    end
end
fclose(fid);

String = strrep(String, '**', '^');
end
