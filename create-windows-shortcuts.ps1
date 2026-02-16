param(
    [string]$ShortcutFolder = [Environment]::GetFolderPath("Desktop")
)

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$shell = New-Object -ComObject WScript.Shell

function New-CmdShortcut {
    param(
        [string]$Name,
        [string]$Target,
        [string]$Arguments = "",
        [string]$WorkDir = $root
    )

    $shortcutPath = Join-Path $ShortcutFolder ("{0}.lnk" -f $Name)
    $shortcut = $shell.CreateShortcut($shortcutPath)
    $shortcut.TargetPath = $Target
    if ($Arguments) {
        $shortcut.Arguments = $Arguments
    }
    $shortcut.WorkingDirectory = $WorkDir
    $shortcut.WindowStyle = 1
    $shortcut.Save()
    Write-Host "Created: $shortcutPath"
}

New-CmdShortcut -Name "UTG-CLI" -Target (Join-Path $root "run-UTG-CLI.bat")
New-CmdShortcut -Name "UTG-WebUI" -Target (Join-Path $root "run-UTG-WebUI.bat")
New-CmdShortcut -Name "DH5aUTG-CLI" -Target (Join-Path $root "run-DH5aUTG-CLI.bat")
New-CmdShortcut -Name "DH5aUTG-WebUI" -Target (Join-Path $root "run-DH5aUTG-WebUI.bat")

Write-Host "All shortcuts are created on: $ShortcutFolder"
