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

# Prevent premature workspace-queue dead-lease expiry during heavy 3D factorizations / swapping (§7.3).
if (-not $env:NYX_QUEUE_LEASE_DEAD_SECS) {
    $env:NYX_QUEUE_LEASE_DEAD_SECS = "60"
}

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

function Build-Workers {
    cargo build --release -p nyx-script
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    cargo build --release -p nyx-engine --bin nyx-worker
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    cargo build --release -p nyx-connect --example fake_mcp
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
            Invoke-Checked "Build nyx-script worker, nyx-worker & nyx-connect fake_mcp" {
                Build-Workers
            }
            Invoke-Checked "nyx-engine contracts, layering, security & vvuq_report gates" {
                cargo test --release -p nyx-engine --test layering --test contracts --test geomdesc --test api_docs --test path_guard --test format --test h5 --test vvuq_phase5_report
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
            Invoke-Checked "Telemetry dry-run smoke check" {
                $py = Get-Command python -ErrorAction SilentlyContinue
                if (-not $py) { $py = Get-Command python3 -ErrorAction SilentlyContinue }
                if ($py) {
                    & $py.Source tools/vvuq_telemetry.py --dry-run
                }
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
            Invoke-Checked "Build nyx-script worker, nyx-worker & nyx-connect fake_mcp" {
                Build-Workers
            }
            Invoke-Checked "Full workspace --release test suite (includes nyx-storage gates, examples & validation)" {
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
            Invoke-Checked "Build nyx-script worker, nyx-worker & nyx-connect fake_mcp" {
                Build-Workers
            }
            Invoke-Checked "Ignored deep verification tests" {
                cargo test --release --workspace -- --ignored
            }
            Invoke-Checked "Windows NYX_THREADS={1,4} & NYX_FACTOR_THREADS={1,4} determinism matrix" {
                foreach ($nt in @("1", "4")) {
                    $env:NYX_THREADS = $nt
                    cargo test --release -p geom --features parallel --test arrangement
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p femcore --test linalg -- determinism
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p femcore --test periodic_mesh3d --test electro_amr --test rt_equilibration --test portable_kernels
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p femcore --features portable-kernels --test portable_kernels
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p nyx-dfs --test qcap_solvers --test farfield_electro
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p nyx-sms --test determinism
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p nyx-ems --test pml_mixed_curvilinear
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p nyx-analysis --test figure3d --test uq_cqed
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    cargo test --release -p nyx-engine --test squadds_engine_smoke --test vvuq_phase4_uq
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    foreach ($crate in @("nyx-dfs", "nyx-ems", "nyx-engine")) {
                        if (Test-Path "crates/$crate/tests/determinism.rs") {
                            cargo test --release -p $crate --test determinism
                            if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                        }
                    }
                    if (Test-Path "crates/nyx-engine/tests/vvuq_phase5_determinism.rs") {
                        cargo test --release -p nyx-engine --test vvuq_phase5_determinism
                        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    }
                }
                $env:NYX_THREADS = "4"
                foreach ($ft in @("1", "4")) {
                    $env:NYX_FACTOR_THREADS = $ft
                    cargo test --release -p nyx-ems --test golden -- golden_cavity_3d_box_v1 golden_bloch_3d_wire_v1 golden_driven_3d_prom_v1
                    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    if (Test-Path "crates/nyx-ems/tests/determinism.rs") {
                        cargo test --release -p nyx-ems --test determinism
                        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    }
                    if (Test-Path "crates/nyx-engine/tests/vvuq_phase5_determinism.rs") {
                        cargo test --release -p nyx-engine --test vvuq_phase5_determinism
                        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                    }
                }
            }
            Invoke-Checked "Full benchmarks/validation/ suites & R22 .fmsh oracle check" {
                cargo test --release -p nyx-engine `
                    --test validation_cqed_cross_code `
                    --test validation_cqed_experiment `
                    --test validation_cqed_tls_scqubits `
                    --test validation_mems_piezo `
                    --test validation_photonics `
                    --test validation_thermal_nafems
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo test --release -p femcore --test contracts -- mesh_fmsh_round_trip_and_version
                if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
                cargo test --release -p nyx-ems --test golden -- golden_soi_wire_and_hybrid_frozen_solver_oracles_r22
            }
            if ($env:VVUQ_SKIP_TELEMETRY -ne "1") {
                $py = Get-Command python -ErrorAction SilentlyContinue
                if (-not $py) { $py = Get-Command python3 -ErrorAction SilentlyContinue }
                if ($py) {
                    $trials = if ($env:VVUQ_TELEMETRY_TRIALS) { $env:VVUQ_TELEMETRY_TRIALS } else { "3" }
                    $outPath = if ($env:VVUQ_TELEMETRY_OUT) { $env:VVUQ_TELEMETRY_OUT } else { "runs/vvuq_telemetry.json" }
                    Invoke-Checked "Quiet-runner best-of-3 telemetry ledger" {
                        & $py.Source tools/vvuq_telemetry.py --trials $trials --out $outPath
                    }
                }
            }
            Write-Host "NOTE: Nightly Miri, Kani (K1-K3), and Bolero LibFuzzer+ASan run on Linux (see tools/vvuq_ci.sh --tier 3)."
        }
    }
    Write-Host "=== [VVUQ Tier $Tier] PASSED ==="
} finally {
    Pop-Location
}
