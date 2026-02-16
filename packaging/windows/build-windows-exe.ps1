param(
    [string]$PythonExe = "",
    [ValidateSet("onefile", "onedir")]
    [string]$Mode = "onefile",
    [switch]$BuildCLI,
    [switch]$BuildInstaller,
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
    param([string]$Command, [string[]]$Args)
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $Command
    $psi.Arguments = ($Args -join " ")
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.WorkingDirectory = $ProjectRoot

    $proc = New-Object System.Diagnostics.Process
    $proc.StartInfo = $psi
    $null = $proc.Start()
    $proc.WaitForExit()
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    Write-Host $stdout
    if ($stderr) { Write-Host $stderr }

    if ($proc.ExitCode -ne 0) {
        throw "Command failed: $Command $($Args -join ' ') (exit $($proc.ExitCode))"
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
    param([string[]]$Args)
    $allArgs = $pythonArgs + $Args
    Invoke-CommandWithOutput -Command $pythonBin -Args $allArgs
}

function Ensure-PyInstaller {
    Write-Host "Installing/refreshing pyinstaller..."
    Invoke-Python @("pip", "install", "--upgrade", "pyinstaller")
}

$oneFileArg = if ($Mode -eq "onefile") { "--onefile" } else { "--onedir" }
$webuiData = '"' + "$ProjectRoot\src;src" + '"'
$webuiBootstrap = '"' + (Join-Path $ProjectRoot "src\webui_bootstrap.py") + '"'
$cliEntry = '"' + (Join-Path $ProjectRoot "src\main.py") + '"'

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

if ($BuildCLI) {
    Write-Host "Build CLI EXE..."
    Invoke-Python @(
        "PyInstaller",
        "--noconfirm",
        "--clean",
        "--name",
        "DH5a-UTG-CLI",
        $oneFileArg,
        "--collect-submodules",
        "src",
        $cliEntry
    )
    $cliBuilt = $true
} elseif (Test-Path "$ProjectRoot\dist\DH5a-UTG-CLI.exe") {
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

Write-Host "Done. Output:"
Get-ChildItem -Path "$ProjectRoot\$OutputDir" -Filter "DH5a-UTG-WebUI*"
