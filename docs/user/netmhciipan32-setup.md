# NetMHCIIpan 3.2 — legacy acquisition and Mac setup

The current standard is 4.3k. This page is specifically for legacy
`predictor_version="3.2"`; see the [current autorun guide](indirect-autorun.md).

## Actual environment and result

This verification host is macOS 26.6.2, arm64. Targeted PATH/common installation
roots and recursive filename searches under Downloads/Desktop found no genuine
3.2 executable or model installation; only a controlled mock was identified.
The search log records pruning and one unrelated unreadable NoMachine log folder.
This is scoped absence, not a claim that every file on the Mac was inspected.
**REAL_PREDICTOR_E2E_NOT_TESTED.** No 4.x substitution or vendor download occurred.

## Official acquisition (the user performs the license step)

1. Open [DTU NetMHCIIpan 3.2](https://services.healthtech.dtu.dk/services/NetMHCIIpan-3.2/).
   Under **Software Downloads**, choose **Version 3.2 → Darwin**; avoid the 4.x entries.
2. The [3.2 Darwin academic request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?packageversion=3.2&platform=Darwin&software=netMHCIIpan&version=3.2)
   requests name, position, institutional email and affiliation. Read the agreement;
   choose acceptance only if entitled and willing. The form emails a download link;
   it also offers license-only review. Noneligible users should contact DTU.
3. The service page separately links [matching Darwin data](https://services.healthtech.dtu.dk/services/NetMHCIIpan-3.2/data.Darwin.tar.gz).
   Obtain matching files under applicable terms and follow the supplied installation
   instructions. Retain original archives and compute SHA-256 locally.
4. The agreement limits use/redistribution and modifications; sections 2 and 7
   should be reviewed before sharing code/data or publishing vendor benchmarks.
   Neoantimon2 does not accept these terms, send the form, redistribute binaries,
   model/training data, or grant a third-party license.

## Apple Silicon compatibility — verify, do not assume

The official 3.2 download selector says **Darwin**, without an arm64 label. The
newer 4.x arm64 listing is not evidence that 3.2 is native arm64. After legitimate
acquisition, inspect the wrapper and its documented native executable with `file`.
This host's existing `arch -x86_64 /usr/bin/true` probe succeeded, but that does not
prove a particular vendor binary/library stack works.

- If DTU supplies a compatible arm64 or universal 3.2 build, use its documented setup.
- For a supported x86_64 Darwin build, verify its documented dependencies and
  execution through the already available Rosetta environment. No Rosetta installation
  or vendor modification is performed by Neoantimon2.
- An i386-only, incompatible or unavailable build needs clarification/a compatible
  **3.2** distribution from DTU at `health-software@dtu.dk`. Do not silently use 4.x.
  Linux is a separate official platform, not a claim that the Mac runner automatically
  manages Docker/virtualization or that a Linux binary executes natively on this Mac.

Suggested local checks after acquisition (replace paths):

```sh
file "/path/to/installed/netMHCIIpan"
"/path/to/installed/netMHCIIpan" -h > netmhciipan32-help.txt 2> netmhciipan32-help.stderr
shasum -a 256 "/path/to/installed/netMHCIIpan"
```

No genuine binary was available to execute these commands here. Retain its original
help and installation documentation before constructing a real invocation.

## Binding to the autorun API

Supply the physical installed model root containing `data/`, an explicit model
release identity, and the path to the legitimate executable. Internal model
symlinks are rejected by this version. The runner checks 3.2-family identity,
`-f`, `-a`, `-inptype`, `-rdir`, supported DR alleles and `-list`/`-listMHC` help;
then runs peptide input in bounded batches. If a vendor wrapper differs, stop
and retain diagnostics rather than modifying its scientific behavior or relabeling
another version. CLI/model compatibility remains unverified until actual execution.

The caller supplies `authorization=list(status="authorized_local_installation",
terms="description of my applicable accepted local terms")`; this records an
explicit declaration, not automated legal verification. Reference setup is the
versioned normalized `hla-reference-1.0.0` manifest with exact literal precursor
proteins and SHA-256. A raw FASTA alone is not that manifest.

Use the [one-command guide](indirect-autorun.md). Raw donor and self inputs,
predictions, conversion, core checks, scoped self exclusion and final exports are
automatic after these prerequisites. The first real technical check should use
the existing small non-patient fixture and retain raw stdout, affinity/core/offset
mapping, model/executable hashes, same-presenter self comparison, count and saved
result reopening. Its count is not required to equal the mock's artificial count.
No official PIRCHE or clinical validation is implied.
