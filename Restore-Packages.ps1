<#
.SYNOPSIS
    Restore NuGet packages for the Amalgam project

.DESCRIPTION
    This script specifically handles NuGet package restoration for C++ projects like Amalgam
    that use packages.config or PackageReference for dependencies like Boost.

.PARAMETER ProjectPath
    Path to the project or solution file. Default is "Amalgam.sln"

.PARAMETER PackagesConfigPath
    Path to packages.config file if different from default location

.PARAMETER Force
    Force reinstall of all packages

.EXAMPLE
    .\Restore-Packages.ps1
    
.EXAMPLE
    .\Restore-Packages.ps1 -ProjectPath "Amalgam\Amalgam.vcxproj" -Force
#>

[CmdletBinding()]
param(
    [string]$ProjectPath = "Amalgam.sln",
    [string]$PackagesConfigPath,
    [switch]$Force
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

# Function to download nuget.exe if not available
function Get-NuGetExecutable {
    $nugetPath = ".\nuget.exe"
    
    if (Test-Path $nugetPath) {
        return $nugetPath
    }
    
    Write-ColorOutput "Downloading nuget.exe..." "Yellow"
    try {
        $nugetUrl = "https://dist.nuget.org/win-x86-commandline/latest/nuget.exe"
        Invoke-WebRequest -Uri $nugetUrl -OutFile $nugetPath
        Write-ColorOutput "NuGet.exe downloaded successfully" "Green"
        return $nugetPath
    } catch {
        Write-ColorOutput "Failed to download nuget.exe: $($_.Exception.Message)" "Red"
        return $null
    }
}

# Function to find packages.config files
function Find-PackagesConfig {
    $configFiles = @()
    
    # Look for packages.config in common locations
    $searchPaths = @(
        "packages.config",
        "Amalgam\packages.config",
        ".\packages.config"
    )
    
    # Also search recursively for packages.config files
    $foundFiles = Get-ChildItem -Path "." -Name "packages.config" -Recurse -ErrorAction SilentlyContinue
    
    foreach ($file in $foundFiles) {
        $configFiles += $file
    }
    
    foreach ($path in $searchPaths) {
        if (Test-Path $path) {
            $configFiles += $path
        }
    }
    
    return $configFiles | Select-Object -Unique
}

try {
    Write-ColorOutput "=== NuGet Package Restoration for Amalgam ===" "Cyan"
    
    # Validate project file exists
    if (-not (Test-Path $ProjectPath)) {
        throw "Project file not found: $ProjectPath"
    }
    
    Write-ColorOutput "Project: $ProjectPath" "Green"
    
    # Find packages.config files
    $packagesConfigs = Find-PackagesConfig
    if ($PackagesConfigPath) {
        $packagesConfigs += $PackagesConfigPath
        $packagesConfigs = $packagesConfigs | Select-Object -Unique
    }
    
    if ($packagesConfigs) {
        Write-ColorOutput "Found packages.config files:" "Green"
        foreach ($config in $packagesConfigs) {
            Write-ColorOutput "  - $config" "Gray"
        }
    }
    
    # Try different restoration methods
    $restoreSuccess = $false
    
    # Method 1: Try nuget.exe restore (best for C++ projects)
    Write-ColorOutput "`nAttempting NuGet package restoration..." "Yellow"
    
    $nugetExe = $null
    if (Test-CommandExists "nuget") {
        $nugetExe = "nuget"
        Write-ColorOutput "Using system nuget.exe" "Gray"
    } else {
        $nugetExe = Get-NuGetExecutable
        if (-not $nugetExe) {
            Write-ColorOutput "Could not obtain nuget.exe" "Red"
        }
    }
    
    if ($nugetExe) {
        Write-ColorOutput "Running: $nugetExe restore $ProjectPath" "Gray"
        
        $restoreArgs = @("restore", $ProjectPath)
        if ($Force) {
            $restoreArgs += "-Force"
        }
        
        & $nugetExe @restoreArgs
        
        if ($LASTEXITCODE -eq 0) {
            $restoreSuccess = $true
            Write-ColorOutput "NuGet restore completed successfully" "Green"
        } else {
            Write-ColorOutput "NuGet restore failed with exit code $LASTEXITCODE" "Red"
        }
    }
    
    # Method 2: Try dotnet restore (for newer project formats)
    if (-not $restoreSuccess -and (Test-CommandExists "dotnet")) {
        Write-ColorOutput "`nTrying dotnet restore..." "Yellow"
        Write-ColorOutput "Running: dotnet restore $ProjectPath" "Gray"
        
        & dotnet restore $ProjectPath
        
        if ($LASTEXITCODE -eq 0) {
            $restoreSuccess = $true
            Write-ColorOutput "dotnet restore completed successfully" "Green"
        } else {
            Write-ColorOutput "dotnet restore failed with exit code $LASTEXITCODE" "Yellow"
        }
    }
    
    # Method 3: Try MSBuild restore target
    if (-not $restoreSuccess) {
        # Try to find MSBuild
        $msbuildPaths = @(
            "${env:ProgramFiles}\Microsoft Visual Studio\2022\Enterprise\MSBuild\Current\Bin\MSBuild.exe",
            "${env:ProgramFiles}\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe",
            "${env:ProgramFiles}\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe",
            "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\Enterprise\MSBuild\Current\Bin\MSBuild.exe",
            "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe",
            "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe"
        )
        
        $msbuildPath = $null
        foreach ($path in $msbuildPaths) {
            if (Test-Path $path) {
                $msbuildPath = $path
                break
            }
        }
        
        if ($msbuildPath) {
            Write-ColorOutput "`nTrying MSBuild restore target..." "Yellow"
            Write-ColorOutput "Running: $msbuildPath $ProjectPath /t:Restore" "Gray"
            
            & $msbuildPath $ProjectPath /t:Restore /v:minimal
            
            if ($LASTEXITCODE -eq 0) {
                $restoreSuccess = $true
                Write-ColorOutput "MSBuild restore completed successfully" "Green"
            } else {
                Write-ColorOutput "MSBuild restore failed with exit code $LASTEXITCODE" "Yellow"
            }
        }
    }
    
    # Check if packages directory was created
    if (Test-Path "packages") {
        Write-ColorOutput "`nPackages directory found" "Green"
        $packageDirs = Get-ChildItem -Path "packages" -Directory | Select-Object -First 5
        Write-ColorOutput "Installed packages:" "Green"
        foreach ($dir in $packageDirs) {
            Write-ColorOutput "  - $($dir.Name)" "Gray"
        }
        if ((Get-ChildItem -Path "packages" -Directory).Count -gt 5) {
            Write-ColorOutput "  - ... and $((Get-ChildItem -Path "packages" -Directory).Count - 5) more" "Gray"
        }
    }
    
    # Specific check for Boost package
    $boostPackage = Get-ChildItem -Path "packages" -Directory -Filter "boost*" -ErrorAction SilentlyContinue
    if ($boostPackage) {
        Write-ColorOutput "Boost package found: $($boostPackage.Name)" "Green"
        
        # Check if boost.targets file exists
        $boostTargets = Get-ChildItem -Path $boostPackage.FullName -Filter "boost.targets" -Recurse -ErrorAction SilentlyContinue
        if ($boostTargets) {
            Write-ColorOutput "boost.targets file found at: $($boostTargets.FullName)" "Green"
        } else {
            Write-ColorOutput "boost.targets file not found in Boost package" "Yellow"
        }
    } else {
        Write-ColorOutput "Boost package not found in packages directory" "Yellow"
    }
    
    if ($restoreSuccess) {
        Write-ColorOutput "`n=== Package Restoration Complete ===" "Green"
        Write-ColorOutput "You can now try building the project with:" "Green"
        Write-ColorOutput "  .\Build.ps1 -ProjectPath `"$ProjectPath`" -Configuration Release -Platform x64" "Gray"
    } else {
        Write-ColorOutput "`n=== Package Restoration Failed ===" "Red"
        Write-ColorOutput "Manual steps to try:" "Yellow"
        Write-ColorOutput "1. Open Visual Studio and right-click solution -> Restore NuGet Packages" "Gray"
        Write-ColorOutput "2. Install NuGet CLI and run: nuget restore $ProjectPath" "Gray"
        Write-ColorOutput "3. Check if packages.config exists and has correct package references" "Gray"
        Write-ColorOutput "4. Verify internet connection and NuGet package sources" "Gray"
        
        # Show packages.config content if found
        if ($packagesConfigs) {
            Write-ColorOutput "`nPackages.config content:" "Yellow"
            foreach ($config in $packagesConfigs) {
                if (Test-Path $config) {
                    Write-ColorOutput "--- $config ---" "Gray"
                    Get-Content $config | Select-Object -First 10 | ForEach-Object {
                        Write-ColorOutput "  $_" "Gray"
                    }
                }
            }
        }
    }
    
} catch {
    Write-ColorOutput "Error: $($_.Exception.Message)" "Red"
    exit 1
}