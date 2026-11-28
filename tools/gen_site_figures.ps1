<#
.SYNOPSIS
    Reproducibly generates the 5 clean showcase figures in website/static/figures/
    via headless nyx scripting (`nyx run ... --workspace ...`).
.EXAMPLE
    .\tools\gen_site_figures.ps1
    .\tools\gen_site_figures.ps1 -NyxExe .\target\release\nyx.exe
#>
[CmdletBinding()]
param(
    [string]$NyxExe = $env:NYX_BIN
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $PSScriptRoot
$ExamplesDir = Join-Path $RepoRoot "examples"
$OutDir = Join-Path $RepoRoot "website\static\figures"

if (-not $NyxExe) {
    $ReleaseNyx = Join-Path $RepoRoot "target\release\nyx.exe"
    if (Test-Path $ReleaseNyx) {
        $NyxExe = $ReleaseNyx
    } else {
        $Cmd = Get-Command nyx -ErrorAction SilentlyContinue
        if ($Cmd) {
            $NyxExe = $Cmd.Source
        } else {
            throw "nyx executable not found. Build with 'cargo build --release -p nyx' or pass -NyxExe."
        }
    }
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$WorkDir = Join-Path ([System.IO.Path]::GetTempPath()) ("nyx-site-figures-" + [System.Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null

try {
    $Examples = @(
        "15_soi_wire.nyx.json",
        "18_wgm_cylindrical_pml.nyx.json",
        "11_anchor_loss.nyx.json",
        "20_cavity_3d_spherical_pml.nyx.json"
    )
    foreach ($ex in $Examples) {
        Copy-Item -LiteralPath (Join-Path $ExamplesDir $ex) -Destination (Join-Path $WorkDir $ex) -Force
    }

    $ScriptPath = Join-Path $WorkDir "render_site_figures.nys"
    $ScriptBody = @'
# Render the 5 showcase PNG plates for website/static/figures/

nyx.project_open(path="15_soi_wire.nyx.json")
nyx.solve_run()
print("15_soi_wire:", nyx.result_info(physics="ems"))
nyx.figure(physics="ems", kind="gallery", columns=2, max_modes=2, scale=2, path="soi_wire_modes.png")

nyx.project_open(path="18_wgm_cylindrical_pml.nyx.json")
nyx.solve_run()
print("18_wgm_cylindrical_pml:", nyx.result_info(physics="ems"))
nyx.figure(physics="ems", kind="gallery", columns=4, max_modes=4, scale=2, path="wgm_cyl_pml.png")

nyx.project_open(path="11_anchor_loss.nyx.json")
nyx.solve_run()
print("11_anchor_loss:", nyx.result_info(physics="sms"))
nyx.figure(physics="sms", kind="gallery", columns=3, max_modes=6, scale=2, path="anchor_loss_modes.png")

nyx.project_open(path="20_cavity_3d_spherical_pml.nyx.json")
nyx.solve_run()
print("20_cavity_3d_spherical_pml:", nyx.result_info(physics="ems"))
nyx.figure(physics="ems", kind="cutaway3d", scale=2, show_edges="regions", path="cavity_3d_cutaway.png")
nyx.figure(physics="ems", kind="slice3d", plane="z", scale=2, show_edges="regions", path="cavity_3d_slice.png")
'@
    [System.IO.File]::WriteAllText($ScriptPath, $ScriptBody, [System.Text.UTF8Encoding]::new($false))

    & $NyxExe run $ScriptPath --workspace $WorkDir --seconds 120
    if ($LASTEXITCODE -ne 0) {
        throw "nyx run exited with code $LASTEXITCODE"
    }

    $Plates = @(
        "soi_wire_modes.png",
        "wgm_cyl_pml.png",
        "anchor_loss_modes.png",
        "cavity_3d_cutaway.png",
        "cavity_3d_slice.png"
    )
    foreach ($p in $Plates) {
        $Src = Join-Path $WorkDir $p
        if (-not (Test-Path $Src)) {
            throw "Expected generated figure missing: $Src"
        }
        Copy-Item -LiteralPath $Src -Destination (Join-Path $OutDir $p) -Force
    }
    Write-Host "Generated 5 showcase plates in $OutDir"
} finally {
    if (Test-Path $WorkDir) {
        Remove-Item -LiteralPath $WorkDir -Recurse -Force -ErrorAction SilentlyContinue
    }
}
