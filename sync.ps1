# Sync this workspace to GitHub: stage, commit and push.
# Usage:
#   .\sync.ps1
#   .\sync.ps1 -Message "fix phase calculation"
param(
    [string]$Message = ""
)

Set-Location $PSScriptRoot

$status = git status --porcelain
if (-not $status) {
    Write-Host "Nothing to sync. Working tree is clean." -ForegroundColor Green
    exit 0
}

git add -A
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

# Safety check: GitHub rejects files over 100 MB and warns above 50 MB.
$tooBig = @()
git diff --cached --name-only | ForEach-Object {
    $path = Join-Path (Get-Location) $_
    if ((Test-Path $path) -and ((Get-Item $path).Length -gt 50MB)) {
        $tooBig += ("{0} ({1:N1} MB)" -f $_, ((Get-Item $path).Length / 1MB))
    }
}

if ($tooBig.Count -gt 0) {
    git reset | Out-Null
    Write-Host "Aborted: these files are larger than 50 MB:" -ForegroundColor Red
    $tooBig | ForEach-Object { Write-Host "  $_" }
    Write-Host "Add them to .gitignore, or use Git LFS." -ForegroundColor Yellow
    exit 1
}

if (-not $Message) {
    $Message = "Sync " + (Get-Date -Format "yyyy-MM-dd HH:mm")
}

git commit -m $Message
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

git push 2>&1 | ForEach-Object { Write-Host $_ }
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "Synced to https://github.com/hellowangxifeng/code" -ForegroundColor Green
