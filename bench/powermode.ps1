# Windows power mode (the "Power mode" overlay of Settings) for a benchmark session, and a sampler of the
# processor's performance (percent of nominal frequency, about once per second).
# The overlay GUIDs are those Windows logged when the mode was changed in Settings (bench/raw/power_mode_events.txt):
#   961cc777-2547-4f9d-8174-7d86181b8a7a 'Best power efficiency' (every table of the paper)
#   00000000-0000-0000-0000-000000000000 'Balanced'
#   ded574b5-45a0-4f42-8737-46345c09c238 'Best performance'
# usage: powermode.ps1 get
#        powermode.ps1 set {efficiency|balanced|performance|GUID}
#        powermode.ps1 sample OUTFILE STOPFILE   (appends one line per sample until STOPFILE exists)
param([string]$Cmd, [string]$A, [string]$B)
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Windows.Forms
Add-Type -Namespace BitdpBench -Name Power -MemberDefinition @'
[DllImport("powrprof.dll")] public static extern uint PowerGetEffectiveOverlayScheme(out Guid g);
[DllImport("powrprof.dll")] public static extern uint PowerSetActiveOverlayScheme(Guid g);
'@
$modes = @{
    efficiency  = '961cc777-2547-4f9d-8174-7d86181b8a7a'
    balanced    = '00000000-0000-0000-0000-000000000000'
    performance = 'ded574b5-45a0-4f42-8737-46345c09c238'
}

function Get-Overlay {
    $g = [Guid]::Empty
    [void][BitdpBench.Power]::PowerGetEffectiveOverlayScheme([ref]$g)
    $name = ($modes.GetEnumerator() | Where-Object { $_.Value -eq "$g" }).Key
    $ac = [System.Windows.Forms.SystemInformation]::PowerStatus.PowerLineStatus -eq 'Online'
    "overlay=$g mode=$name ac=$ac"
}

switch ($Cmd) {
    'get' { Get-Overlay }
    'set' {
        $g = if ($modes.ContainsKey($A)) { $modes[$A] } else { $A }
        $rc = [BitdpBench.Power]::PowerSetActiveOverlayScheme([Guid]$g)
        if ($rc -ne 0) { throw "PowerSetActiveOverlayScheme returned $rc" }
        Get-Overlay
    }
    'sample' {
        # Get-Counter waits one second for the rate counter, which sets the sampling period.
        while (-not (Test-Path $B)) {
            $v = (Get-Counter '\Processor Information(_Total)\% Processor Performance').CounterSamples[0].CookedValue
            "$(Get-Date -Format o) perf=$([math]::Round($v, 1)) $(Get-Overlay)" | Add-Content -Path $A
        }
    }
    default { throw 'usage: powermode.ps1 get | set MODE | sample OUTFILE STOPFILE' }
}
