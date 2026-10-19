param(
    [Parameter(Mandatory)][string]$Archive,
    [string]$Checksum = "$Archive.sha256"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Assert($ok, [string]$message) {
    if (!$ok) { throw $message }
}

$archivePath = (Resolve-Path -LiteralPath $Archive).Path
$checksumPath = (Resolve-Path -LiteralPath $Checksum).Path
$line = (Get-Content -LiteralPath $checksumPath -Raw).Trim()
Assert ($line -match '^(?<hash>[0-9A-Fa-f]{64})\s{2}(?<name>[^\\/]+)$') `
    'Archive checksum must be HASH, two spaces, and the archive filename.'
Assert ($matches.name -ceq [IO.Path]::GetFileName($archivePath)) `
    'Checksum filename does not match the archive.'
Assert ((Get-FileHash -LiteralPath $archivePath -Algorithm SHA256).Hash -ceq
        $matches.hash.ToUpperInvariant()) 'Archive SHA-256 mismatch.'

$tempRoot = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
$temp = Join-Path $tempRoot ("sparse_poreg_verify_" + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $temp | Out-Null
try {
    Expand-Archive -LiteralPath $archivePath -DestinationPath $temp
    $manifest = Join-Path $temp 'SHA256SUMS'
    Assert (Test-Path -LiteralPath $manifest -PathType Leaf) 'Archive lacks SHA256SUMS.'
    $seen = @{}
    foreach ($entry in Get-Content -LiteralPath $manifest) {
        Assert ($entry -match '^(?<hash>[0-9A-Fa-f]{64})\s{2}(?<path>.+)$') `
            "Malformed SHA256SUMS entry: $entry"
        $wantHash = $matches.hash.ToUpperInvariant()
        $rel = $matches.path.Replace('/', [IO.Path]::DirectorySeparatorChar)
        Assert (![IO.Path]::IsPathRooted($rel) -and $rel -notmatch '(^|[\\/])\.\.([\\/]|$)') `
            "Unsafe SHA256SUMS path: $rel"
        $path = [IO.Path]::GetFullPath((Join-Path $temp $rel))
        Assert ($path.StartsWith($temp + [IO.Path]::DirectorySeparatorChar,
                [StringComparison]::OrdinalIgnoreCase)) "Unsafe SHA256SUMS path: $rel"
        Assert (!$seen.ContainsKey($path)) "Duplicate SHA256SUMS path: $rel"
        $seen[$path] = 1
        Assert (Test-Path -LiteralPath $path -PathType Leaf) "Archive file is missing: $rel"
        Assert ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash -ceq
                $wantHash) "Archive member SHA-256 mismatch: $rel"
    }
    $actual = @(Get-ChildItem -LiteralPath $temp -Recurse -File |
        Where-Object { $_.FullName -cne $manifest }).Count
    Assert ($seen.Count -eq $actual) 'SHA256SUMS does not cover every archive member.'
}
finally {
    $resolved = [IO.Path]::GetFullPath($temp)
    if ($resolved.StartsWith($tempRoot, [StringComparison]::OrdinalIgnoreCase) -and
            [IO.Path]::GetFileName($resolved) -match '^sparse_poreg_verify_[0-9a-f]{32}$') {
        Remove-Item -LiteralPath $resolved -Recurse -Force -ErrorAction SilentlyContinue
    }
}

Write-Host "PASS: $([IO.Path]::GetFileName($archivePath)) and all archived files match SHA-256."
