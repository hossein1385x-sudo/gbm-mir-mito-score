<#
  publish_to_github.ps1 - one-time helper to push this repository to GitHub.
  Edit the three variables below, then run from PowerShell:
      powershell -ExecutionPolicy Bypass -File .\publish_to_github.ps1
  No token or password is stored here: Git asks for your GitHub login in the browser
  (Git Credential Manager) the first time you push.
#>

# ---------- EDIT THESE ----------
$RepoDir    = "C:\Users\pishtazrayne- HEMATI\Desktop\rapo\gbm-mir-mito-score\gbm-mir-mito-score"
$GitHubUser = "hossein1385x-sudo"
$RepoName   = "gbm-mir-mito-score"
$GitName    = "Seyed Hossein Hosseini"
$GitEmail   = "hossein1385x@gmail.com"
# --------------------------------

$ErrorActionPreference = "Stop"

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw "Git is not installed. Download it from https://git-scm.com/download/win"
}
if (-not (Test-Path (Join-Path $RepoDir "README.md"))) {
    throw "README.md not found in $RepoDir - check `$RepoDir."
}
Set-Location $RepoDir

# 1) Safety check: refuse to continue if data or large files are present in the folder
$bad = Get-ChildItem -Recurse -File -Force |
       Where-Object { $_.FullName -notmatch '\\\.git\\' -and
                      ( $_.Extension -in ".rds",".RData",".rda",".docx",".gz",".tsv",".csv",".xlsx" -or
                        $_.Length -gt 5MB ) }
if ($bad) {
    Write-Host "These files should NOT go to GitHub (data / generated / large):" -ForegroundColor Red
    $bad | ForEach-Object { Write-Host "  $($_.FullName)  ($([math]::Round($_.Length/1MB,2)) MB)" }
    throw "Remove them (or add them to .gitignore) and run again."
}
if (Test-Path ".Renviron") { Write-Host ".Renviron exists locally; it is git-ignored and will not be pushed." }

# 2) Initialise and commit
if (-not (Test-Path ".git")) { git init }
git config user.name  $GitName
git config user.email $GitEmail
git branch -M main
git add .
Write-Host "`n=== Files that will be committed ===" -ForegroundColor Cyan
git status --short
$answer = Read-Host "`nLooks right? Type YES to commit and push"
if ($answer -ne "YES") { Write-Host "Stopped. Nothing was pushed."; exit }

git commit -m "Initial release: miR-Mito Score pipeline (TCGA-GBM, CGGA validation)"

# 3) Push (create the EMPTY repo on github.com first, without README/license)
$remote = "https://github.com/$GitHubUser/$RepoName.git"
if (-not (git remote)) { git remote add origin $remote } else { git remote set-url origin $remote }
git push -u origin main

Write-Host "`nDone: https://github.com/$GitHubUser/$RepoName" -ForegroundColor Green
