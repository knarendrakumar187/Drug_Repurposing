# Run one or more Paper 2 task scripts with a single Neo4j password prompt.
# The password is kept only in this PowerShell process for the duration of the run.
#   & .\Paper2_Results\code\run_tasks.ps1 task3_druggability task4_drug_coverage
param([Parameter(Mandatory = $true, ValueFromRemainingArguments = $true)][string[]]$Tasks)

$code = $PSScriptRoot
if (-not $env:NEO4J_PASSWORD) {
    $secure = Read-Host -AsSecureString "Neo4j password"
    $env:NEO4J_PASSWORD = [Runtime.InteropServices.Marshal]::PtrToStringAuto(
        [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secure))
}
$env:PYTHONIOENCODING = "utf-8"
try {
    foreach ($t in $Tasks) {
        Write-Host "===== $t ====="
        python (Join-Path $code "$t.py")
        if ($LASTEXITCODE -ne 0) { Write-Host "$t failed (exit $LASTEXITCODE); stopping."; break }
        Write-Host "===== $t finished ====="
    }
}
finally {
    Remove-Item Env:NEO4J_PASSWORD -ErrorAction SilentlyContinue
}
