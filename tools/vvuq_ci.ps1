<#
.SYNOPSIS
    VVUQ 3-tier verification runner (Windows PowerShell mirror of .github/workflows/vvuq-ci.yml).
.EXAMPLE
    .\tools\vvuq_ci.ps1 -Tier 1
    .\tools\vvuq_ci.ps1 -Tier 2
    .\tools\vvuq_ci.ps1 -Tier 3
#>
[CmdletBinding()]
param(
    [ValidateSet(1, 2, 3)]
    [int]$Tier = 1,
    [switch]$ExcludeNyx
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $RepoRoot

function Invoke-Checked {
    param(
        [string]$Label,
        [scriptblock]$Command
    )
    Write-Host "=== [VVUQ Tier $Tier] $Label ==="
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Step failed ($LASTEXITCODE): $Label"
    }
}

try {
    switch ($Tier) {
        1 {
            Invoke-Checked "ADR 65 formatting (7 crates)" {
                cargo fmt --check -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine
            }
            if ($ExcludeNyx) {
                Invoke-Checked "Clippy (--release, --exclude nyx)" {
                    cargo clippy --release --workspace --exclude nyx --all-targets --no-deps -- -D warnings
                }
            } else {
                Invoke-Checked "Clippy (--release, workspace)" {
                    cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
                }
            }
            Invoke-Checked "Build nyx-script worker & nyx-connect fake_mcp" {
                cargo build --release -p nyx-script
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo build --release -p nyx-connect --example fake_mcp
            }
            Invoke-Checked "nyx-engine contracts, layering & security gates" {
                cargo test --release -p nyx-engine --test layering --test contracts --test geomdesc --test api_docs --test path_guard --test format --test h5
            }
            Invoke-Checked "femcore contracts" {
                cargo test --release -p femcore --test contracts --test linalg_contract --test material_contract
            }
            Invoke-Checked "Workspace crate contracts (9 crates)" {
                cargo test --release -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-store -p nyx-mcp -p nyx-connect -p nyx-install --test contracts
            }
            Invoke-Checked "nyx-script contract" {
                cargo test --release -p nyx-script --test contract
            }
            Invoke-Checked "P0/P1 soundness, SSRF & proptest regressions" {
                cargo test --release -p nyx-engine --lib material_fetch
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo test --release -p geom --test props --test gds
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo test --release -p nyx-sms --test finite_guard
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                foreach ($crate in @("nyx-dfs", "nyx-ems")) {
                    if (Test-Path "crates/$crate/tests/finite_guard.rs") {
                        cargo test --release -p $crate --test finite_guard
                        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    }
                }
                foreach ($target in @("fuzz_gds", "fuzz_nyx_storage", "fuzz_script_trace", "fuzz_path_and_archive", "fuzz_mcp_framing")) {
                    $hit = Get-ChildItem -Path "crates" -Filter "$target.rs" -Recurse -Depth 3 -ErrorAction SilentlyContinue | Select-Object -First 1
                    if ($hit) {
                        $pkg = Split-Path -Leaf (Split-Path -Parent (Split-Path -Parent $hit.FullName))
                        cargo test --release -p $pkg --test $target
                        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    }
                }
            }
        }
        2 {
            Invoke-Checked "Build nyx-script worker & nyx-connect fake_mcp" {
                cargo build --release -p nyx-script
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo build --release -p nyx-connect --example fake_mcp
            }
            Invoke-Checked "Full workspace --release test suite (includes nyx-storage gates & examples)" {
                cargo test --release --workspace
            }
            Invoke-Checked "Differential portable-kernels (default features)" {
                cargo test --release -p femcore --test portable_kernels
            }
            Invoke-Checked "Differential portable-kernels (--features portable-kernels)" {
                cargo test --release -p femcore --features portable-kernels --test portable_kernels
            }
            Invoke-Checked "femcore linalg_contract (--features portable-kernels)" {
                cargo test --release -p femcore --features portable-kernels --test linalg_contract
            }
        }
        3 {
            Invoke-Checked "Build nyx-script worker & nyx-connect fake_mcp" {
                cargo build --release -p nyx-script
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo build --release -p nyx-connect --example fake_mcp
            }
            Invoke-Checked "Ignored deep verification tests" {
                cargo test --release --workspace -- --ignored
            }
            Invoke-Checked "Windows determinism & portable-kernels parity" {
                $env:NYX_THREADS = "1"
                cargo test --release -p nyx-sms --test determinism
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo test --release -p femcore --test portable_kernels
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                $env:NYX_THREADS = "4"
                cargo test --release -p nyx-sms --test determinism
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo test --release -p femcore --features portable-kernels --test portable_kernels
            }
            Write-Host "NOTE: Nightly Miri, Kani (K1-K3), and Bolero LibFuzzer+ASan run on Linux (see tools/vvuq_ci.sh --tier 3)."
        }
    }
    Write-Host "=== [VVUQ Tier $Tier] PASSED ==="
} finally {
    Pop-Location
}
