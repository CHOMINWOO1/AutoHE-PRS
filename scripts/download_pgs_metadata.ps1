$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
$OutDir = Join-Path $Root "data\pgs_catalog\metadata"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$Files = @(
    "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/pgs_scores_list.txt",
    "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_scores.csv"
)

foreach ($Url in $Files) {
    $Name = Split-Path -Leaf $Url
    $Target = Join-Path $OutDir $Name
    if (Test-Path -LiteralPath $Target) {
        Write-Host "Already exists: $Target"
        continue
    }
    Write-Host "Downloading $Url"
    & curl.exe --ssl-no-revoke -L --fail --retry 3 --output $Target $Url
    if ($LASTEXITCODE -ne 0) {
        throw "curl.exe failed with exit code $LASTEXITCODE"
    }
}

Write-Host "Downloaded metadata to $OutDir"

