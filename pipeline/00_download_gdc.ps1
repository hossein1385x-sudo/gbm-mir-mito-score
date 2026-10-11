<#
.SYNOPSIS
  Download one part of a GDC manifest (RNA-seq or miRNA-seq) with gdc-client, resumably.

.DESCRIPTION
  The GDC manifest was split into numbered parts (rnaseq_part_N.tsv / mirna_part_N.tsv)
  inside WorkDir. This script:
    * skips files that are already complete (final file present and no .partial file),
    * writes a temporary manifest containing only the pending files,
    * runs gdc-client with 4 parallel connections (so one throttled file does not block the queue),
    * re-checks how many files are complete afterwards.
  Re-running the script resumes where it stopped; it never starts from scratch.

.PARAMETER DataType
  "rnaseq" or "mirna".

.PARAMETER Part
  Number of the manifest part to download.

.PARAMETER WorkDir
  Folder that contains the manifest parts. Downloads go to <WorkDir>\<DataType>\<file uuid>\.
  Use the same folder as GBM_WORK_DIR (see config.R).

.PARAMETER GdcClient
  Path to gdc-client.exe (default: gdc-client.exe found on PATH).

.EXAMPLE
  .\pipeline\00_download_gdc.ps1 -DataType rnaseq -Part 8 -WorkDir D:\gbm_work
  .\pipeline\00_download_gdc.ps1 -DataType mirna  -Part 6 -WorkDir D:\gbm_work -GdcClient C:\tools\gdc-client.exe
#>
param(
    [Parameter(Mandatory = $true)][ValidateSet("rnaseq", "mirna")][string]$DataType,
    [Parameter(Mandatory = $true)][int]$Part,
    [string]$WorkDir   = $(if ($env:GBM_WORK_DIR) { $env:GBM_WORK_DIR } else { (Join-Path (Get-Location) "work") }),
    [string]$GdcClient = "gdc-client.exe"
)

$ErrorActionPreference = "Stop"
$dest = Join-Path $WorkDir $DataType
New-Item -ItemType Directory -Force -Path $dest | Out-Null
Set-Location $WorkDir

function Get-PendingManifestLines {
    param([string]$ManifestPath, [string]$DestPath)
    $allLines   = Get-Content -Path $ManifestPath -Encoding UTF8
    $header     = $allLines[0]
    $dataLines  = $allLines[1..($allLines.Count - 1)]
    $headerCols = $header -split "`t"
    $idIndex    = [array]::IndexOf($headerCols, "id")
    $nameIndex  = [array]::IndexOf($headerCols, "filename")
    if ($idIndex -eq -1 -or $nameIndex -eq -1) { throw "Columns 'id' or 'filename' not found in the manifest." }
    $pending = New-Object System.Collections.Generic.List[string]
    $skipped = 0
    foreach ($line in $dataLines) {
        if ([string]::IsNullOrWhiteSpace($line)) { continue }
        $cols        = $line -split "`t"
        $finalFile   = Join-Path (Join-Path $DestPath $cols[$idIndex]) $cols[$nameIndex]
        $partialFile = "$finalFile.partial"
        # Only complete files (final file present AND no .partial file) are skipped.
        if ((Test-Path $finalFile) -and -not (Test-Path $partialFile)) { $skipped++ }
        else { $pending.Add($line) }
    }
    return [PSCustomObject]@{ Header = $header; Pending = $pending; Skipped = $skipped; Total = $dataLines.Count }
}

function Write-ManifestNoBom {
    param([string]$Path, [string[]]$Lines)
    $utf8NoBom = New-Object System.Text.UTF8Encoding $false
    [System.IO.File]::WriteAllText($Path, (($Lines -join "`n") + "`n"), $utf8NoBom)
}

$manifestPath = Join-Path $WorkDir "${DataType}_part_$Part.tsv"
$result = Get-PendingManifestLines -ManifestPath $manifestPath -DestPath $dest
Write-Host "Total: $($result.Total) | already complete: $($result.Skipped) | remaining: $($result.Pending.Count)"

if ($result.Pending.Count -eq 0) {
    Write-Host "All files of this part are already downloaded. Nothing to do." -ForegroundColor Green
    return
}

$tempManifest = Join-Path ([System.IO.Path]::GetTempPath()) "pending_${DataType}_part_$Part.tsv"
Write-ManifestNoBom -Path $tempManifest -Lines (@($result.Header) + $result.Pending)

Write-Host "`n--- Temporary manifest (first 3 lines) ---" -ForegroundColor Cyan
Get-Content $tempManifest -TotalCount 3

Write-Host "`n--- Running gdc-client with 4 parallel connections ---" -ForegroundColor Yellow
& $GdcClient download -d $dest -m $tempManifest -n 4 --retry-amount 5 --wait-time 3 --http-chunk-size 1048576
Write-Host "`ngdc-client exit code: $LASTEXITCODE" -ForegroundColor Magenta

$after = Get-PendingManifestLines -ManifestPath $manifestPath -DestPath $dest
Write-Host "`nAfter this run: $($after.Skipped)/$($after.Total) complete, $($after.Pending.Count) remaining." -ForegroundColor Cyan
if ($after.Pending.Count -gt 0) {
    Write-Host "Run the same command again to continue; it resumes, it does not restart." -ForegroundColor Yellow
}
