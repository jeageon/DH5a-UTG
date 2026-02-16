#define MyAppName "DH5a-UTG"
#define MyAppVersion "1.0.0"
#define MyOutputDir ".\\installer"

[Setup]
AppId={{A7E8F8F0-5B0A-4A6C-B8CF-1A4A0E2B2A11}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
OutputDir={#MyOutputDir}
OutputBaseFilename=DH5a-UTG-Setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern
DisableProgramGroupPage=yes
PrivilegesRequired=lowest

[Languages]
Name: "korean"; MessagesFile: "compiler:Languages\\Korean.isl"

[Files]
Source: "..\\..\\dist\\DH5a-UTG-WebUI.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\\..\\dist\\DH5a-UTG-CLI.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\\..\\README.md"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\\DH5a-UTG-WebUI"; Filename: "{app}\\DH5a-UTG-WebUI.exe"
Name: "{group}\\DH5a-UTG-CLI"; Filename: "{app}\\DH5a-UTG-CLI.exe"
Name: "{autodesktop}\\DH5a-UTG-WebUI"; Filename: "{app}\\DH5a-UTG-WebUI.exe"
Name: "{autodesktop}\\DH5a-UTG-CLI"; Filename: "{app}\\DH5a-UTG-CLI.exe"

[Run]
Filename: "{app}\\DH5a-UTG-WebUI.exe"; Description: "DH5a-UTG WebUI 실행"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{app}"

