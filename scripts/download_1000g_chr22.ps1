$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
$OutDir = Join-Path $Root "data\1000genomes\phase3"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$Files = @(
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
)

foreach ($Url in $Files) {
    $Name = Split-Path -Leaf $Url
    $Target = Join-Path $OutDir $Name
    if (Test-Path -LiteralPath $Target) {
        Write-Host "Already exists: $Target"
        continue
    }
    Write-Host "Downloading $Name"
    Invoke-WebRequest -Uri $Url -OutFile $Target
}

Write-Host "Downloaded files to $OutDir"

