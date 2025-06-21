<#
.SYNOPSIS
    PowerShell build script for MSBuild projects outside of Visual Studio

.DESCRIPTION
    This script provides a comprehensive solution for building MSBuild projects from the command line
    without requiring Visual Studio. It supports various build configurations, platforms, and options
    commonly used in MSBuild projects.

.PARAMETER ProjectPath
    Path to the project file (.csproj, .vcxproj, .sln, etc.) to build

.PARAMETER Configuration
    Build configuration (Debug, Release, etc.). Default is "Release"

.PARAMETER Platform
    Target platform (x86, x64, AnyCPU, etc.). Default is "AnyCPU"

.PARAMETER OutputPath
    Custom output directory for build artifacts

.PARAMETER Clean
    Clean the project before building

.PARAMETER Rebuild
    Rebuild the project (clean + build)

.PARAMETER Restore
    Restore NuGet packages before building

.PARAMETER Verbosity
    MSBuild verbosity level (quiet, minimal, normal, detailed, diagnostic). Default is "normal"

.PARAMETER Properties
    Additional MSBuild properties as hashtable

.PARAMETER Targets
    Specific targets to build (comma-separated)

.PARAMETER MaxCpuCount
    Maximum number of concurrent processes to use during build

.PARAMETER BinaryLog
    Enable binary logging and specify log file path

.PARAMETER FileLogger
    Enable file logging and specify log file path

.PARAMETER NoLogo
    Suppress MSBuild startup banner

.PARAMETER WarnAsError
    Treat warnings as errors

.PARAMETER UseDotNet
    Use 'dotnet build' instead of MSBuild.exe

.EXAMPLE
    .\Build.ps1 -ProjectPath "MyProject.csproj"
    
.EXAMPLE
    .\Build.ps1 -ProjectPath "MySolution.sln" -Configuration "Debug" -Platform "x64" -Clean

.EXAMPLE
    .\Build.ps1 -ProjectPath "MyProject.csproj" -UseDotNet -BinaryLog "build.binlog"

.EXAMPLE
    .\Build.ps1 -ProjectPath "MyProject.csproj" -Properties @{DefineConstants="CUSTOM_BUILD";WarningLevel=4}
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectPath,
    
    [string]$Configuration = "Release",
    
    [string]$Platform = "AnyCPU",
    
    [string]$OutputPath,
    
    [switch]$Clean,
    
    [switch]$Rebuild,
    
    [switch]$Restore,
    
    [ValidateSet("quiet", "minimal", "normal", "detailed", "diagnostic")]
    [string]$Verbosity = "normal",
    
    [hashtable]$Properties = @{},
    
    [string[]]$Targets,
    
    [int]$MaxCpuCount,
    
    [string]$BinaryLog,
    
    [string]$FileLogger,
    
    [switch]$NoLogo,
    
    [switch]$WarnAsError,
    
    [switch]$UseDotNet
)

# Function to write colored output
function Write-ColorOutput {
    param(
        [string]$Message,
        [string]$Color = "White"
    )
    Write-Host $Message -ForegroundColor $Color
}

# Function to check if a command exists
function Test-CommandExists {
    param([string]$Command)
    $null = Get-Command $Command -ErrorAction SilentlyContinue
    return $?
}

# Function to find MSBuild executable
function Find-MSBuild {
    # Try to find MSBuild in common locations
    $msbuildPaths = @(
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\Enterprise\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\2019\Enterprise\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\2019\Professional\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\2019\Community\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\Enterprise\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2019\Enterprise\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2019\Professional\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2019\Community\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Shared\MSBuild\Current\Bin\MSBuild.exe",
        "${env:ProgramFiles}\dotnet\sdk\*\MSBuild.dll"
    )
    
    foreach ($path in $msbuildPaths) {
        if (Test-Path $path) {
            return $path
        }
    }
    
    # Try to find MSBuild in PATH
    if (Test-CommandExists "msbuild") {
        return "msbuild"
    }
    
    return $null
}

# Main script execution
try {
    Write-ColorOutput "=== MSBuild PowerShell Build Script ===" "Cyan"
    Write-ColorOutput "Project: $ProjectPath" "Green"
    Write-ColorOutput "Configuration: $Configuration" "Green"
    Write-ColorOutput "Platform: $Platform" "Green"
    
    # Validate project file exists
    if (-not (Test-Path $ProjectPath)) {
        throw "Project file not found: $ProjectPath"
    }
    
    # Determine build tool to use
    $buildCommand = ""
    $buildArgs = @()
    
    if ($UseDotNet) {
        # Use dotnet CLI
        if (-not (Test-CommandExists "dotnet")) {
            throw ".NET CLI (dotnet) not found. Please install .NET SDK or use -UseDotNet:$false"
        }
        
        Write-ColorOutput "Using .NET CLI (dotnet build)" "Yellow"
        $buildCommand = "dotnet"
        
        # Determine dotnet command based on operation
        if ($Clean -and -not $Rebuild) {
            $buildArgs += "clean"
        } elseif ($Rebuild) {
            # For rebuild, we'll clean first then build
            $buildArgs += "build"
            $buildArgs += "--no-incremental"
        } else {
            $buildArgs += "build"
        }
        
        $buildArgs += $ProjectPath
        
        # Add configuration
        $buildArgs += "--configuration", $Configuration
        
        # Add verbosity
        $buildArgs += "--verbosity", $Verbosity
        
        # Add restore if requested
        if ($Restore) {
            $buildArgs += "--no-restore:false"
        } else {
            $buildArgs += "--no-restore"
        }
        
        # Add no-logo if requested
        if ($NoLogo) {
            $buildArgs += "--nologo"
        }
        
        # Add binary log
        if ($BinaryLog) {
            $buildArgs += "-bl:$BinaryLog"
        }
        
        # Add properties
        $allProperties = $Properties.Clone()
        if ($Platform -and $Platform -ne "AnyCPU") {
            $allProperties["Platform"] = $Platform
        }
        if ($OutputPath) {
            $allProperties["OutputPath"] = $OutputPath
        }
        if ($WarnAsError) {
            $allProperties["TreatWarningsAsErrors"] = "true"
        }
        
        foreach ($prop in $allProperties.GetEnumerator()) {
            $buildArgs += "-p:$($prop.Key)=$($prop.Value)"
        }
        
    } else {
        # Use MSBuild.exe
        $msbuildPath = Find-MSBuild
        if (-not $msbuildPath) {
            throw "MSBuild not found. Please install Visual Studio Build Tools or use -UseDotNet switch"
        }
        
        Write-ColorOutput "Using MSBuild: $msbuildPath" "Yellow"
        $buildCommand = $msbuildPath
        
        $buildArgs += $ProjectPath
        
        # Add targets
        if ($Targets) {
            $buildArgs += "/t:$($Targets -join ';')"
        } elseif ($Clean -and -not $Rebuild) {
            $buildArgs += "/t:Clean"
        } elseif ($Rebuild) {
            $buildArgs += "/t:Rebuild"
        } else {
            $buildArgs += "/t:Build"
        }
        
        # Add properties
        $allProperties = $Properties.Clone()
        $allProperties["Configuration"] = $Configuration
        if ($Platform -and $Platform -ne "AnyCPU") {
            $allProperties["Platform"] = $Platform
        }
        if ($OutputPath) {
            $allProperties["OutputPath"] = $OutputPath
        }
        if ($WarnAsError) {
            $allProperties["TreatWarningsAsErrors"] = "true"
        }
        
        foreach ($prop in $allProperties.GetEnumerator()) {
            $buildArgs += "/p:$($prop.Key)=$($prop.Value)"
        }
        
        # Add verbosity
        $buildArgs += "/v:$Verbosity"
        
        # Add max CPU count
        if ($MaxCpuCount) {
            $buildArgs += "/m:$MaxCpuCount"
        } else {
            $buildArgs += "/m"
        }
        
        # Add no-logo if requested
        if ($NoLogo) {
            $buildArgs += "/nologo"
        }
        
        # Add binary log
        if ($BinaryLog) {
            $buildArgs += "/bl:$BinaryLog"
        }
        
        # Add file logger
        if ($FileLogger) {
            $buildArgs += "/fl", "/flp:logfile=$FileLogger"
        }
    }
    
    # Restore packages if requested and not using dotnet (dotnet handles restore automatically)
    if ($Restore -and -not $UseDotNet) {
        Write-ColorOutput "Restoring NuGet packages..." "Yellow"
        
        # Try different NuGet restore methods
        $restoreSuccess = $false
        
        # Method 1: Try nuget.exe restore (best for C++ projects with packages.config)
        if (Test-CommandExists "nuget") {
            Write-ColorOutput "Using nuget.exe for package restore..." "Gray"
            & nuget restore $ProjectPath
            if ($LASTEXITCODE -eq 0) {
                $restoreSuccess = $true
            }
        }
        
        # Method 2: Try dotnet restore (for PackageReference style projects)
        if (-not $restoreSuccess -and (Test-CommandExists "dotnet")) {
            Write-ColorOutput "Using dotnet restore for package restore..." "Gray"
            & dotnet restore $ProjectPath
            if ($LASTEXITCODE -eq 0) {
                $restoreSuccess = $true
            }
        }
        
        # Method 3: Try MSBuild restore target
        if (-not $restoreSuccess -and $msbuildPath) {
            Write-ColorOutput "Using MSBuild restore target..." "Gray"
            & $msbuildPath $ProjectPath /t:Restore /v:minimal
            if ($LASTEXITCODE -eq 0) {
                $restoreSuccess = $true
            }
        }
        
        if (-not $restoreSuccess) {
            Write-ColorOutput "Warning: Package restore failed or no restore tool found" "Yellow"
            Write-ColorOutput "You may need to manually restore packages using:" "Yellow"
            Write-ColorOutput "  - nuget restore $ProjectPath" "Gray"
            Write-ColorOutput "  - dotnet restore $ProjectPath" "Gray"
            Write-ColorOutput "  - Or restore packages in Visual Studio" "Gray"
        }
    }
    
    # Clean if requested and doing rebuild with dotnet
    if ($Rebuild -and $UseDotNet) {
        Write-ColorOutput "Cleaning project..." "Yellow"
        & dotnet clean $ProjectPath --configuration $Configuration --verbosity $Verbosity
        if ($LASTEXITCODE -ne 0) {
            throw "Clean failed"
        }
    }
    
    # Execute build command
    Write-ColorOutput "Executing build command..." "Yellow"
    Write-ColorOutput "Command: $buildCommand $($buildArgs -join ' ')" "Gray"
    
    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
    
    & $buildCommand @buildArgs
    $exitCode = $LASTEXITCODE
    
    $stopwatch.Stop()
    $buildTime = $stopwatch.Elapsed
    
    # Report results
    if ($exitCode -eq 0) {
        Write-ColorOutput "Build succeeded in $($buildTime.ToString('mm\:ss\.fff'))" "Green"
        
        # Display output information
        if ($OutputPath) {
            Write-ColorOutput "Output directory: $OutputPath" "Green"
        }
        
        if ($BinaryLog) {
            Write-ColorOutput "Binary log saved to: $BinaryLog" "Green"
        }
        
        if ($FileLogger) {
            Write-ColorOutput "Build log saved to: $FileLogger" "Green"
        }
        
    } else {
        Write-ColorOutput "Build failed with exit code $exitCode after $($buildTime.ToString('mm\:ss\.fff'))" "Red"
        exit $exitCode
    }
    
} catch {
    Write-ColorOutput "Error: $($_.Exception.Message)" "Red"
    exit 1
}

# Additional helper functions for advanced scenarios

<#
.SYNOPSIS
    Helper function to build multiple projects in sequence
#>
function Build-MultipleProjects {
    param(
        [string[]]$ProjectPaths,
        [hashtable]$BuildParameters = @{}
    )
    
    foreach ($project in $ProjectPaths) {
        Write-ColorOutput "Building project: $project" "Cyan"
        $params = $BuildParameters.Clone()
        $params["ProjectPath"] = $project
        
        & $PSCommandPath @params
        if ($LASTEXITCODE -ne 0) {
            throw "Build failed for project: $project"
        }
    }
}

<#
.SYNOPSIS
    Helper function to setup build environment variables
#>
function Set-BuildEnvironment {
    param(
        [hashtable]$EnvironmentVariables = @{}
    )
    
    # Common MSBuild environment variables
    $defaultEnvVars = @{
        "MSBUILDTERMINALLOGGER" = "auto"
        "MSBUILDLOGIMPORTS" = "1"
        "MSBUILDLOGALLTASKEVENTS" = "1"
    }
    
    # Merge with user-provided variables
    $allEnvVars = $defaultEnvVars + $EnvironmentVariables
    
    foreach ($envVar in $allEnvVars.GetEnumerator()) {
        [Environment]::SetEnvironmentVariable($envVar.Key, $envVar.Value, "Process")
        Write-ColorOutput "Set environment variable: $($envVar.Key) = $($envVar.Value)" "Gray"
    }
}
