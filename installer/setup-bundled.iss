; Script generated by the Inno Script Studio Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "DTOcean Hydrodynamic Data"
#define MyAppNameSafe "dtocean-hydrodynamic-data"
#define MyAppVersion "1.0.0"
#define MyAppPublisher "DTOcean"
#define MyAppURL "http://www.dtocean.eu/"
#define BinPath "bin"
#define WaveIncludePath "data\dtocean_wec"
#define TidalIncludePath "data\dtocean_tidal"
#define IniDirName "DTOcean Hydrodynamics"

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{54C4421C-EB34-4D00-A9F3-B37593213A75}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
;AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DisableReadyPage=yes
DefaultDirName={sd}\{#MyAppPublisher}
DisableDirPage=yes
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
DisableFinishedPage=yes
OutputBaseFilename={#MyAppNameSafe}-{#MyAppVersion}
Compression=lzma
SolidCompression=yes
PrivilegesRequired=lowest
UninstallFilesDir={commonappdata}\{#MyAppPublisher}\{#MyAppName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: "..\{#BinPath}\*"; DestDir: "{code:GetPath}\{#BinPath}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\{#WaveIncludePath}\*"; DestDir: "{code:GetPath}\{#WaveIncludePath}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\{#TidalIncludePath}\*"; DestDir: "{code:GetPath}\{#TidalIncludePath}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\dtocean_hydro\config\install.ini"; DestDir: "{commonappdata}\{#MyAppPublisher}\{#IniDirName}"; DestName: "install.ini"; Permissions: users-modify
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"

[INI]
Filename: "{commonappdata}\{#MyAppPublisher}\{#IniDirName}\install.ini"; Section: "dtocean_tidal"; Key: "include_path"; String: "{code:GetPath}\{#TidalIncludePath}";
Filename: "{commonappdata}\{#MyAppPublisher}\{#IniDirName}\install.ini"; Section: "dtocean_wec"; Key: "bin_path"; String: "{code:GetPath}\{#BinPath}";
Filename: "{commonappdata}\{#MyAppPublisher}\{#IniDirName}\install.ini"; Section: "dtocean_wec"; Key: "include_path"; String: "{code:GetPath}\{#WaveIncludePath}";

[Code]
function GetPath(Value: String): String;
var
  OrigPath: string;
begin
  if RegQueryStringValue(HKCU, 'Software\DTOcean', 'INSTALL_DIR', OrigPath) then
    Result := OrigPath
  else
    MsgBox('DTOcean installation directory not found!', mbError, MB_OK);
end;

const
  BN_CLICKED = 0;
  WM_COMMAND = $0111;
  CN_BASE = $BC00;
  CN_COMMAND = CN_BASE + WM_COMMAND;

procedure CurPageChanged(CurPageID: Integer);
var
  Param: Longint;
begin
  { if we are on the ready page, then... }
  if CurPageID = wpReady then
  begin
    { the result of this is 0, just to be precise... }
    Param := 0 or BN_CLICKED shl 16;
    { post the click notification message to the next button }
    PostMessage(WizardForm.NextButton.Handle, CN_COMMAND, Param, 0);
  end;
end;