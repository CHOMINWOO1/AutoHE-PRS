param(
    [Parameter(Mandatory = $true)]
    [string]$PgsId,

    [ValidateSet("formatted", "GRCh37", "GRCh38")]
    [string]$Build = "GRCh37"
)

$ErrorActionPreference = "Stop"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$Root = Split-Path -Parent $PSScriptRoot
$OutDir = Join-Path $Root "data\pgs_catalog\$PgsId"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

if ($Build -eq "formatted") {
    $FileName = "$PgsId.txt.gz"
    $Url = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/$PgsId/ScoringFiles/$FileName"
} else {
    $FileName = "${PgsId}_hmPOS_${Build}.txt.gz"
    $Url = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/$PgsId/ScoringFiles/Harmonized/$FileName"
}

$Target = Join-Path $OutDir $FileName
if (Test-Path -LiteralPath $Target) {
    Write-Host "Already exists: $Target"
} else {
    Write-Host "Downloading $Url"
    try {
        Invoke-WebRequest -Uri $Url -OutFile $Target
    } catch {
        Write-Host "Invoke-WebRequest failed; retrying with curl.exe"
        & curl.exe --ssl-no-revoke -L --fail --retry 3 --output $Target $Url
        if ($LASTEXITCODE -ne 0) {
            throw "curl.exe failed with exit code $LASTEXITCODE"
        }
    }
}

Write-Host "Saved to $Target"
