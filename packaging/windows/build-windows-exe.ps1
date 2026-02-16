param(
    [string]$PythonExe = "",
    [ValidateSet("onefile", "onedir")]
    [string]$Mode = "onefile",
    [switch]$BuildCLI,
    [switch]$BuildInstaller,
    [switch]$CreateDesktopShortcuts = $true,
    [string]$ShortcutFolder = "",
    [string]$OutputDir = "dist",
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
)

$ErrorActionPreference = "Stop"

function Resolve-PythonExecutable {
    param([string]$Preferred)
    if ($Preferred -and (Test-Path $Preferred)) {
        return $Preferred
    }

    $venvPython = Join-Path $ProjectRoot ".venv\Scripts\python.exe"
    if (Test-Path $venvPython) {
        return $venvPython
    }

    $pyPath = Get-Command py -ErrorAction SilentlyContinue
    if ($pyPath) {
        return "py -3"
    }

    $sysPython = Get-Command python -ErrorAction SilentlyContinue
    if ($sysPython) {
        return $sysPython.Source
    }

    throw "Python을 찾지 못했습니다. Python 3 설치 또는 .venv를 확인하세요."
}

function Invoke-CommandWithOutput {
    param([string]$Command, [string[]]$CommandArgs)
    Write-Host "Running: $Command $($CommandArgs -join ' ')"

    $oldLocation = Get-Location
    try {
        Set-Location $ProjectRoot
        & $Command @CommandArgs 2>&1 | ForEach-Object { Write-Host $_ }
    }
    finally {
        Set-Location $oldLocation
    }

    if ($LASTEXITCODE -ne 0) {
        throw "Command failed: $Command $($CommandArgs -join ' ') (exit $LASTEXITCODE)"
    }
}

$pythonExe = Resolve-PythonExecutable -Preferred $PythonExe
Write-Host "Using Python: $pythonExe"

if ($pythonExe -eq "py -3") {
    $pythonBin = "py"
    $pythonArgs = @("-3", "-m")
} else {
    $pythonBin = $pythonExe
    $pythonArgs = @("-m")
}

function Invoke-Python {
    param([string[]]$PyArgs)
    $allArgs = $pythonArgs + $PyArgs
    Invoke-CommandWithOutput -Command $pythonBin -CommandArgs $allArgs
}

function Ensure-PyInstaller {
    Write-Host "Installing/refreshing pyinstaller..."
    Invoke-Python @("pip", "install", "--upgrade", "pyinstaller")
}

$oneFileArg = if ($Mode -eq "onefile") { "--onefile" } else { "--onedir" }
$webuiData = "$ProjectRoot\src;src"
$webuiBootstrap = Join-Path $ProjectRoot "src\webui_bootstrap.py"
$cliEntry = Join-Path $ProjectRoot "src\main.py"
$desktopShortcutFolder = if ($ShortcutFolder) { $ShortcutFolder } else { [Environment]::GetFolderPath("Desktop") }

function New-DesktopShortcut {
    param(
        [string]$Name,
        [string]$Target,
        [string]$WorkDir = $ProjectRoot
    )

    if (-not (Test-Path $Target)) {
        throw "Target executable not found: $Target"
    }

    $shell = New-Object -ComObject WScript.Shell
    $shortcutPath = Join-Path $desktopShortcutFolder ("$Name.lnk")
    $shortcut = $shell.CreateShortcut($shortcutPath)
    $shortcut.TargetPath = $Target
    $shortcut.WorkingDirectory = $WorkDir
    $shortcut.WindowStyle = 1
    $shortcut.Description = "DH5a-UTG"
    $shortcut.Save()
    Write-Host "Created desktop shortcut: $shortcutPath"
}

function New-ShortcutSet {
    param(
        [bool]$BuildCLIExecutable
    )

    $webuiExe = Resolve-ExePath -Name "DH5a-UTG-WebUI"
    New-DesktopShortcut -Name "DH5aUTG-WebUI" -Target $webuiExe -WorkDir $ProjectRoot

    if ($BuildCLIExecutable) {
        $cliExe = Resolve-ExePath -Name "DH5a-UTG-CLI"
        New-DesktopShortcut -Name "DH5aUTG-CLI" -Target $cliExe -WorkDir $ProjectRoot
    }
}

function Resolve-ExePath {
    param(
        [Parameter(Mandatory=$true)]
        [string]$Name
    )

    $expected = @(
        (Join-Path $ProjectRoot (Join-Path $OutputDir ("$Name.exe"))),
        (Join-Path $ProjectRoot (Join-Path $OutputDir (Join-Path $Name "$Name.exe")))
    )

    foreach ($candidate in $expected) {
        if (Test-Path $candidate) {
            return $candidate
        }
    }

    Write-Host "Expected path not found for $Name. Scanning workspace for fallback match..."
    $fallback = Get-ChildItem -Path $ProjectRoot -Recurse -File -Filter "$Name*.exe" -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending

    if (-not $fallback) {
        Write-Host "No candidate found with pattern '$Name*.exe'."
        if (Test-Path (Join-Path $ProjectRoot $OutputDir)) {
            Write-Host "Dist folder contents:"
            Get-ChildItem -Path (Join-Path $ProjectRoot $OutputDir) -Recurse -File | ForEach-Object { Write-Host $_.FullName }
        } else {
            Write-Host "Dist folder does not exist: $(Join-Path $ProjectRoot $OutputDir)"
        }
        return ""
    }

    return $fallback[0].FullName
}

function Assert-Executable {
    param(
        [Parameter(Mandatory=$true)]
        [string]$Name
    )

    $path = Resolve-ExePath -Name $Name
    if (-not $path) {
        throw "Missing .\\$OutputDir\\$Name.exe (or .\\$OutputDir\\$Name\\$Name.exe)"
    }
    return $path
}

Ensure-PyInstaller

if (-not (Test-Path "$ProjectRoot\$OutputDir")) {
    New-Item -ItemType Directory -Path "$ProjectRoot\$OutputDir" | Out-Null
}

Write-Host "Build WebUI EXE..."
$cliBuilt = $false
Invoke-Python @(
    "PyInstaller",
    "--noconfirm",
    "--clean",
    "--distpath",
    (Join-Path $ProjectRoot $OutputDir),
    "--name",
    "DH5a-UTG-WebUI",
    $oneFileArg,
    "--collect-submodules",
    "src",
    "--collect-all",
    "streamlit",
    "--add-data",
    $webuiData,
    $webuiBootstrap
)
$webuiExe = Assert-Executable -Name "DH5a-UTG-WebUI"

if ($BuildCLI) {
    Write-Host "Build CLI EXE..."
    Invoke-Python @(
        "PyInstaller",
        "--noconfirm",
        "--clean",
        "--distpath",
        (Join-Path $ProjectRoot $OutputDir),
        "--name",
        "DH5a-UTG-CLI",
        $oneFileArg,
        "--collect-submodules",
        "src",
        $cliEntry
    )
    $cliExe = Assert-Executable -Name "DH5a-UTG-CLI"
    $cliBuilt = $true
} elseif (Test-Path (Join-Path $ProjectRoot (Join-Path $OutputDir "DH5a-UTG-CLI.exe"))) {
    $cliBuilt = $true
}

if ($BuildInstaller) {
    if (-not $cliBuilt) {
        throw "BuildInstaller requires CLI executable. Use -BuildCLI for this release."
    }

    $iscc = Get-Command iscc -ErrorAction SilentlyContinue
    if (-not $iscc) {
        Write-Warning "iscc(시스템 PATH)가 없어 Inno Setup 설치 파일을 만들 수 없습니다."
        Write-Host "Inno Setup을 설치한 뒤 다음으로 반복 실행하세요: -BuildInstaller"
        Write-Host "Install script: packaging\\windows\\dh5a_utg_setup.iss"
    } else {
        & iscc (Join-Path $PSScriptRoot "dh5a_utg_setup.iss")
    }
}

if ($CreateDesktopShortcuts) {
    New-ShortcutSet -BuildCLIExecutable $cliBuilt
}

Write-Host "Workspace exe scan (DH5a-UTG-*.exe):"
Get-ChildItem -Path $ProjectRoot -Recurse -File -Filter "DH5a-UTG-*.exe" -ErrorAction SilentlyContinue |
    Sort-Object FullName |
    ForEach-Object { Write-Host $_.FullName }

Write-Host "Done. Output:"
Get-ChildItem -Path (Join-Path $ProjectRoot $OutputDir) -Recurse | Where-Object { $_.Name -like "*DH5a-UTG-*.exe" } | ForEach-Object { Write-Host $_.FullName }
