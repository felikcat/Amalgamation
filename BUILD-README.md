# MSBuild PowerShell Build Scripts

This directory contains PowerShell scripts for building MSBuild projects outside of Visual Studio, specifically designed for the Amalgam project but adaptable to other MSBuild-based projects.

## Files

- **[`Build.ps1`](Build.ps1)** - Main PowerShell build script with comprehensive MSBuild support
- **[`Restore-Packages.ps1`](Restore-Packages.ps1)** - NuGet package restoration script for C++ projects
- **[`Build-Amalgam-Example.ps1`](Build-Amalgam-Example.ps1)** - Example script showing various build scenarios for the Amalgam project
- **[`BUILD-README.md`](BUILD-README.md)** - This documentation file

## Prerequisites

### Required
- **PowerShell 5.1 or later** (Windows PowerShell or PowerShell Core)
- **MSBuild** - One of the following:
  - Visual Studio 2019/2022 (any edition)
  - Visual Studio Build Tools 2019/2022
  - .NET SDK (for .NET projects)

### Optional
- **.NET CLI** (`dotnet`) - For .NET projects and package restoration
- **Git** - For version control operations

## Quick Start

### Basic Usage

```powershell
# IMPORTANT: Restore NuGet packages first (required for C++ projects with NuGet dependencies)
.\Restore-Packages.ps1 -ProjectPath "Amalgam.sln"

# Build the Amalgam solution in Release configuration
.\Build.ps1 -ProjectPath "Amalgam.sln" -Configuration "Release" -Platform "x64" -Restore

# Clean and rebuild in Debug configuration
.\Build.ps1 -ProjectPath "Amalgam.sln" -Configuration "Debug" -Platform "x64" -Rebuild -Restore

# Build with binary logging for troubleshooting
.\Build.ps1 -ProjectPath "Amalgam.sln" -BinaryLog "build.binlog" -Restore
```

### Run Examples

```powershell
# Run all example build scenarios (includes package restoration)
.\Build-Amalgam-Example.ps1
```

## Script Parameters

### [`Build.ps1`](Build.ps1) Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ProjectPath` | String | *Required* | Path to project file (.sln, .csproj, .vcxproj, etc.) |
| `Configuration` | String | "Release" | Build configuration (Debug, Release, etc.) |
| `Platform` | String | "AnyCPU" | Target platform (x86, x64, AnyCPU, etc.) |
| `OutputPath` | String | | Custom output directory for build artifacts |
| `Clean` | Switch | | Clean the project before building |
| `Rebuild` | Switch | | Rebuild the project (clean + build) |
| `Restore` | Switch | | Restore NuGet packages before building |
| `Verbosity` | String | "normal" | MSBuild verbosity (quiet, minimal, normal, detailed, diagnostic) |
| `Properties` | Hashtable | @{} | Additional MSBuild properties |
| `Targets` | String[] | | Specific targets to build |
| `MaxCpuCount` | Int | | Maximum concurrent processes |
| `BinaryLog` | String | | Enable binary logging to specified file |
| `FileLogger` | String | | Enable file logging to specified file |
| `NoLogo` | Switch | | Suppress MSBuild startup banner |
| `WarnAsError` | Switch | | Treat warnings as errors |
| `UseDotNet` | Switch | | Use `dotnet build` instead of MSBuild.exe |

## Usage Examples

### 1. Basic Builds

```powershell
# Simple release build
.\Build.ps1 -ProjectPath "Amalgam.sln"

# Debug build with specific platform
.\Build.ps1 -ProjectPath "Amalgam.sln" -Configuration "Debug" -Platform "x64"

# Clean before building
.\Build.ps1 -ProjectPath "Amalgam.sln" -Clean
```

### 2. Advanced Builds

```powershell
# Rebuild with custom properties
.\Build.ps1 -ProjectPath "Amalgam.sln" -Rebuild -Properties @{
    "WarningLevel" = "4"
    "TreatWarningsAsErrors" = "false"
    "MultiProcessorCompilation" = "true"
}

# Build with maximum parallelism
.\Build.ps1 -ProjectPath "Amalgam.sln" -MaxCpuCount 8

# Build specific targets
.\Build.ps1 -ProjectPath "Amalgam.sln" -Targets @("Clean", "Build", "Publish")
```

### 3. Logging and Diagnostics

```powershell
# Binary logging for MSBuild Structured Log Viewer
.\Build.ps1 -ProjectPath "Amalgam.sln" -BinaryLog "build.binlog"

# Text file logging
.\Build.ps1 -ProjectPath "Amalgam.sln" -FileLogger "build.log"

# Detailed verbosity for troubleshooting
.\Build.ps1 -ProjectPath "Amalgam.sln" -Verbosity "detailed"

# Diagnostic verbosity (very verbose)
.\Build.ps1 -ProjectPath "Amalgam.sln" -Verbosity "diagnostic"
```

### 4. .NET Projects

```powershell
# Use dotnet CLI instead of MSBuild.exe
.\Build.ps1 -ProjectPath "MyProject.csproj" -UseDotNet

# Restore packages and build
.\Build.ps1 -ProjectPath "MyProject.csproj" -UseDotNet -Restore
```

## Build Environment Detection

The script automatically detects and uses the best available build tools:

1. **MSBuild Detection Order:**
   - Visual Studio 2022 (Enterprise → Professional → Community)
   - Visual Studio 2019 (Enterprise → Professional → Community)
   - MSBuild in PATH
   - .NET SDK MSBuild (when using `-UseDotNet`)

2. **Platform Detection:**
   - Automatically handles x86/x64 MSBuild selection
   - Supports cross-compilation scenarios

## Troubleshooting

### Common Issues

1. **"MSBuild not found"**
   ```powershell
   # Install Visual Studio Build Tools or use .NET CLI
   .\Build.ps1 -ProjectPath "project.sln" -UseDotNet
   ```

2. **"This project references NuGet package(s) that are missing" (Boost, etc.)**
   ```powershell
   # Use the dedicated package restoration script
   .\Restore-Packages.ps1 -ProjectPath "Amalgam.sln"
   
   # Or restore during build
   .\Build.ps1 -ProjectPath "Amalgam.sln" -Restore
   
   # Force reinstall packages if needed
   .\Restore-Packages.ps1 -ProjectPath "Amalgam.sln" -Force
   ```

3. **Build fails with missing dependencies**
   ```powershell
   # Restore packages first
   .\Build.ps1 -ProjectPath "project.sln" -Restore
   ```

4. **Need detailed error information**
   ```powershell
   # Use binary logging for analysis
   .\Build.ps1 -ProjectPath "project.sln" -BinaryLog "debug.binlog" -Verbosity "diagnostic"
   ```

### Analyzing Binary Logs

Binary logs (`.binlog` files) can be analyzed using:
- **MSBuild Structured Log Viewer** - Download from GitHub
- **Visual Studio** - File → Open → Open Binary Log
- **Command line** - `msbuild build.binlog` to replay the build

### Environment Variables

The script respects common MSBuild environment variables:
- `MSBUILDTERMINALLOGGER` - Controls terminal logger behavior
- `MSBUILDLOGIMPORTS` - Logs import information
- `MSBUILDLOGALLTASKEVENTS` - Logs all task events

## Integration with CI/CD

### Azure DevOps

```yaml
- task: PowerShell@2
  displayName: 'Build Amalgam'
  inputs:
    filePath: 'Build.ps1'
    arguments: '-ProjectPath "Amalgam.sln" -Configuration "$(BuildConfiguration)" -Platform "$(BuildPlatform)" -BinaryLog "$(Agent.TempDirectory)/build.binlog"'
```

### GitHub Actions

```yaml
- name: Build Amalgam
  run: |
    .\Build.ps1 -ProjectPath "Amalgam.sln" -Configuration "Release" -Platform "x64" -BinaryLog "build.binlog"
  shell: pwsh
```

### Jenkins

```groovy
powershell '''
    .\\Build.ps1 -ProjectPath "Amalgam.sln" -Configuration "${env.BUILD_CONFIGURATION}" -Platform "x64" -BinaryLog "build.binlog"
'''
```

## Advanced Features

### Custom Build Functions

The script exports helper functions for advanced scenarios:

```powershell
# Build multiple projects in sequence
Build-MultipleProjects -ProjectPaths @("Project1.csproj", "Project2.csproj") -BuildParameters @{Configuration="Release"}

# Set custom build environment
Set-BuildEnvironment -EnvironmentVariables @{
    "MSBUILDTERMINALLOGGER" = "auto"
    "CUSTOM_BUILD_FLAG" = "true"
}
```

### Extending the Script

You can extend the script by:
1. Adding custom properties to the `$Properties` hashtable
2. Implementing custom pre/post-build steps
3. Adding project-specific logic based on file extensions

## Performance Tips

1. **Use parallel builds:** `-MaxCpuCount` (or let MSBuild auto-detect)
2. **Enable incremental builds:** Avoid `-Clean` unless necessary
3. **Use binary logging:** Minimal performance impact, maximum diagnostics
4. **Optimize verbosity:** Use "minimal" for faster builds, "detailed" for debugging

## Support

For issues specific to:
- **MSBuild:** Check the [MSBuild documentation](https://docs.microsoft.com/en-us/visualstudio/msbuild/)
- **Amalgam project:** Refer to the project's specific documentation
- **This script:** Check the inline help with `Get-Help .\Build.ps1 -Full`

## License

This build script is provided as-is for use with the Amalgam project and similar MSBuild-based projects.