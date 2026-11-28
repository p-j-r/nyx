// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Streaming `NYXPK001` archive reader/writer, PE Authenticode-aware 96-byte `NYXINST1` self-extracting
//! trailer, `.zip` and directory/`Nyx.app` payload sources, and integrity verification.

use crate::manifest::{FileBase, FileRole, TargetOs};
use crate::notices;
use crate::safety::validate_rel_path;
use anyhow::{bail, Context, Result};
use flate2::read::ZlibDecoder;
use flate2::write::ZlibEncoder;
use flate2::Compression;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

pub const PACK_MAGIC: &[u8; 8] = b"NYXPK001";
pub const TRAILER_MAGIC: &[u8; 8] = b"NYXINST1";
pub const TRAILER_VERSION: u8 = 1;
pub const CODEC_ZLIB_RS: u8 = 1;
pub const TRAILER_LEN: usize = 96;

/// Hard decompression limits to neutralize corrupted or malicious archives.
pub const MAX_SINGLE_FILE_BYTES: u64 = 160 * 1024 * 1024;
pub const MAX_TOTAL_PAYLOAD_BYTES: u64 = 512 * 1024 * 1024;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BundleEntry {
    pub rel_path: String,
    pub base: FileBase,
    pub role: FileRole,
    pub uncompressed_size: u64,
    pub compressed_size: u64,
    /// Offset in bytes from the end of the TOC JSON (i.e. start of the compressed data section).
    pub data_offset: u64,
    pub blake3: String,
    pub mode: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BundleToc {
    pub format: u32,
    pub nyx_version: String,
    pub target_triple: String,
    pub os: TargetOs,
    pub entries: Vec<BundleEntry>,
}

/// Fixed-size 96-byte little-endian trailer (padded on an 8-byte boundary before Authenticode certs).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Trailer {
    pub version: u8,
    pub codec: u8,
    pub flags: u16,
    pub toc_len: u32,
    pub payload_offset: u64,
    pub payload_len: u64,
    pub toc_blake3: [u8; 32],
    pub payload_blake3: [u8; 32],
}

impl Trailer {
    pub fn to_bytes(&self) -> [u8; TRAILER_LEN] {
        let mut out = [0u8; TRAILER_LEN];
        out[0..8].copy_from_slice(TRAILER_MAGIC);
        out[8] = self.version;
        out[9] = self.codec;
        out[10..12].copy_from_slice(&self.flags.to_le_bytes());
        out[12..16].copy_from_slice(&self.toc_len.to_le_bytes());
        out[16..24].copy_from_slice(&self.payload_offset.to_le_bytes());
        out[24..32].copy_from_slice(&self.payload_len.to_le_bytes());
        out[32..64].copy_from_slice(&self.toc_blake3);
        out[64..96].copy_from_slice(&self.payload_blake3);
        out
    }

    pub fn from_bytes(buf: &[u8; TRAILER_LEN]) -> Result<Self> {
        if &buf[0..8] != TRAILER_MAGIC {
            bail!("missing NYXINST1 trailer magic");
        }
        let version = buf[8];
        if version != TRAILER_VERSION {
            bail!("unsupported installer trailer version {version}");
        }
        let codec = buf[9];
        if codec != CODEC_ZLIB_RS {
            bail!("unsupported installer payload codec {codec}");
        }
        let flags = u16::from_le_bytes(buf[10..12].try_into().unwrap());
        let toc_len = u32::from_le_bytes(buf[12..16].try_into().unwrap());
        let payload_offset = u64::from_le_bytes(buf[16..24].try_into().unwrap());
        let payload_len = u64::from_le_bytes(buf[24..32].try_into().unwrap());
        let mut toc_blake3 = [0u8; 32];
        toc_blake3.copy_from_slice(&buf[32..64]);
        let mut payload_blake3 = [0u8; 32];
        payload_blake3.copy_from_slice(&buf[64..96]);
        Ok(Self {
            version,
            codec,
            flags,
            toc_len,
            payload_offset,
            payload_len,
            toc_blake3,
            payload_blake3,
        })
    }
}

/// Extension point for release authenticity verification (Phases 1-2 verify BLAKE3 integrity; Phase 3 plugs in Ed25519/Minisign).
pub trait PayloadVerifier {
    fn verify_trailer(&self, trailer: &Trailer, toc_bytes: &[u8]) -> Result<()>;
    fn verify_payload_stream(
        &self,
        trailer: &Trailer,
        actual_payload_blake3: &[u8; 32],
    ) -> Result<()>;
}

#[derive(Default, Clone, Copy, Debug)]
pub struct IntegrityVerifier;

impl PayloadVerifier for IntegrityVerifier {
    fn verify_trailer(&self, trailer: &Trailer, toc_bytes: &[u8]) -> Result<()> {
        let actual = blake3::hash(toc_bytes);
        if actual.as_bytes() != &trailer.toc_blake3 {
            bail!("installer TOC BLAKE3 verification failed (corrupted executable or bundle)");
        }
        Ok(())
    }

    fn verify_payload_stream(
        &self,
        trailer: &Trailer,
        actual_payload_blake3: &[u8; 32],
    ) -> Result<()> {
        if actual_payload_blake3 != &trailer.payload_blake3 {
            bail!("installer payload BLAKE3 verification failed (download appears corrupted)");
        }
        Ok(())
    }
}

/// Inspects a PE header for `DataDirectory[4]` (Authenticode `SECURITY` table).
/// Returns `Some(cert_file_offset)` when a PE certificate table is present at `cert_offset >= 96`.
pub fn pe_security_cert_offset<R: Read + Seek>(r: &mut R) -> Result<Option<u64>> {
    r.seek(SeekFrom::Start(0))?;
    let mut dos = [0u8; 64];
    if r.read_exact(&mut dos).is_err() || &dos[0..2] != b"MZ" {
        return Ok(None);
    }
    let e_lfanew = u32::from_le_bytes(dos[60..64].try_into().unwrap()) as u64;
    if !(64..=1024 * 1024).contains(&e_lfanew) {
        return Ok(None);
    }
    r.seek(SeekFrom::Start(e_lfanew))?;
    let mut pe = [0u8; 24];
    if r.read_exact(&mut pe).is_err() || &pe[0..4] != b"PE\0\0" {
        return Ok(None);
    }
    let opt_start = e_lfanew + 24;
    r.seek(SeekFrom::Start(opt_start))?;
    let mut magic_buf = [0u8; 2];
    if r.read_exact(&mut magic_buf).is_err() {
        return Ok(None);
    }
    let magic = u16::from_le_bytes(magic_buf);
    let (rva_count_off, dd_start) = match magic {
        0x10B => (opt_start + 92, opt_start + 96),
        0x20B => (opt_start + 108, opt_start + 112),
        _ => return Ok(None),
    };
    r.seek(SeekFrom::Start(rva_count_off))?;
    let mut cnt_buf = [0u8; 4];
    r.read_exact(&mut cnt_buf)?;
    let num_rvas = u32::from_le_bytes(cnt_buf);
    if num_rvas <= 4 {
        return Ok(None);
    }
    r.seek(SeekFrom::Start(dd_start + 4 * 8))?;
    let mut sec_dd = [0u8; 8];
    r.read_exact(&mut sec_dd)?;
    let cert_offset = u32::from_le_bytes(sec_dd[0..4].try_into().unwrap());
    let cert_size = u32::from_le_bytes(sec_dd[4..8].try_into().unwrap());
    if cert_size > 0 && (cert_offset as usize) >= TRAILER_LEN {
        Ok(Some(cert_offset as u64))
    } else {
        Ok(None)
    }
}

/// Locates an embedded `NYXINST1` trailer inside an installer executable (Authenticode-signed or unsigned).
/// Returns `Ok(None)` when the executable is a bare stub without an appended payload.
pub fn locate_embedded_trailer<R: Read + Seek>(r: &mut R) -> Result<Option<(Trailer, u64)>> {
    let file_len = r.seek(SeekFrom::End(0))?;
    if file_len < TRAILER_LEN as u64 {
        return Ok(None);
    }
    let trailer_end = match pe_security_cert_offset(r)? {
        Some(off) if off <= file_len => off,
        _ => file_len,
    };
    if trailer_end < TRAILER_LEN as u64 {
        return Ok(None);
    }
    let trailer_start = trailer_end - TRAILER_LEN as u64;
    r.seek(SeekFrom::Start(trailer_start))?;
    let mut buf = [0u8; TRAILER_LEN];
    r.read_exact(&mut buf)?;
    if &buf[0..8] != TRAILER_MAGIC {
        return Ok(None);
    }
    let trailer = Trailer::from_bytes(&buf)?;
    if trailer.payload_offset.saturating_add(trailer.payload_len) > trailer_start {
        bail!("corrupted NYXINST1 trailer offsets");
    }
    Ok(Some((trailer, trailer_start)))
}

#[derive(Clone, Debug)]
pub struct StagedItem {
    pub source_path: PathBuf,
    pub rel_path: String,
    pub base: FileBase,
    pub role: FileRole,
    pub mode: u32,
}

#[derive(Clone, Debug)]
pub struct PackSpec {
    pub nyx_version: String,
    pub target_triple: String,
    pub os: TargetOs,
    pub items: Vec<StagedItem>,
}

pub fn binary_names(include_worker: bool) -> Vec<&'static str> {
    let mut out = vec!["nyx", "nyx-mcp"];
    if include_worker {
        out.push("nyx-worker");
    }
    out.extend(["nyx-script-worker", "nyx-connect"]);
    out
}

pub fn exe_filename(stem: &str, os: TargetOs) -> String {
    if os == TargetOs::Windows {
        format!("{stem}.exe")
    } else {
        stem.to_string()
    }
}

/// Inputs controlling how `collect_pack_spec` stages binaries, assets, and third-party notices.
#[derive(Debug, Clone, Copy)]
pub struct PackCollectOpts<'a> {
    pub bin_source_dir: &'a Path,
    pub repo_root: &'a Path,
    pub uninstaller_stub: &'a Path,
    pub notices_text: Option<&'a str>,
    pub os: TargetOs,
    pub target_triple: &'a str,
    pub nyx_version: &'a str,
    pub include_worker: bool,
}

/// Collects all files from a build binary directory (`CARGO_TARGET_DIR/release`) + repo assets into a deterministic `PackSpec`.
pub fn collect_pack_spec(opts: &PackCollectOpts<'_>) -> Result<(PackSpec, Option<PathBuf>)> {
    let mut items = Vec::new();

    for stem in binary_names(opts.include_worker) {
        let fname = exe_filename(stem, opts.os);
        let src = opts.bin_source_dir.join(&fname);
        if !src.is_file() {
            bail!("missing release binary '{}'", src.display());
        }
        items.push(StagedItem {
            source_path: src,
            rel_path: fname,
            base: FileBase::Bin,
            role: FileRole::Binary,
            mode: 0o755,
        });
    }

    if !opts.uninstaller_stub.is_file() {
        bail!(
            "missing uninstaller stub binary '{}'",
            opts.uninstaller_stub.display()
        );
    }
    items.push(StagedItem {
        source_path: opts.uninstaller_stub.to_path_buf(),
        rel_path: exe_filename("nyx-uninstall", opts.os),
        base: FileBase::Bin,
        role: FileRole::Uninstaller,
        mode: 0o755,
    });

    for (folder, prefix, ext_filter) in [
        (opts.repo_root.join("materials"), "materials", Some("toml")),
        (
            opts.repo_root.join("examples"),
            "examples",
            Some("nyx.json"),
        ),
    ] {
        if !folder.is_dir() {
            bail!("missing asset directory '{}'", folder.display());
        }
        let mut entries: Vec<PathBuf> = std::fs::read_dir(&folder)?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| p.is_file())
            .filter(|p| match ext_filter {
                Some("toml") => p.extension().and_then(|e| e.to_str()) == Some("toml"),
                Some("nyx.json") => p
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.ends_with(".nyx.json")),
                _ => p.file_name().and_then(|n| n.to_str()) != Some("README.md"),
            })
            .collect();
        entries.sort();
        for p in entries {
            let name = p
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or_default()
                .to_string();
            items.push(StagedItem {
                source_path: p,
                rel_path: format!("{prefix}/{name}"),
                base: FileBase::Resource,
                role: FileRole::Asset,
                mode: 0o644,
            });
        }
    }

    for (src, rel) in [
        (opts.repo_root.join("LICENSE"), "LICENSE"),
        (opts.repo_root.join("NOTICE"), "NOTICE"),
    ] {
        if !src.is_file() {
            bail!("missing license file '{}'", src.display());
        }
        items.push(StagedItem {
            source_path: src,
            rel_path: rel.to_string(),
            base: FileBase::Resource,
            role: FileRole::License,
            mode: 0o644,
        });
    }

    let generated_text = match opts.notices_text {
        Some(t) => t.to_string(),
        None => notices::generate_third_party_notices(opts.repo_root, Some(opts.target_triple))?,
    };
    let tmp_notice = std::env::temp_dir().join(format!("nyx-notices-{}.txt", std::process::id()));
    std::fs::write(&tmp_notice, generated_text.as_bytes())?;
    items.push(StagedItem {
        source_path: tmp_notice.clone(),
        rel_path: "THIRD_PARTY_NOTICES.txt".to_string(),
        base: FileBase::Resource,
        role: FileRole::License,
        mode: 0o644,
    });

    // Deterministic order across builds.
    items.sort_by(|a, b| {
        a.base
            .as_cmp_key()
            .cmp(&b.base.as_cmp_key())
            .then(a.rel_path.cmp(&b.rel_path))
    });
    Ok((
        PackSpec {
            nyx_version: opts.nyx_version.to_string(),
            target_triple: opts.target_triple.to_string(),
            os: opts.os,
            items,
        },
        Some(tmp_notice),
    ))
}

/// Output hashes and TOC metadata emitted by `write_nyxpack` for self-extracting trailers.
#[derive(Debug, Clone)]
pub struct WrittenPack {
    pub toc_bytes: Vec<u8>,
    pub payload_len: u64,
    pub toc_blake3: [u8; 32],
    pub payload_blake3: [u8; 32],
}

impl FileBase {
    fn as_cmp_key(self) -> u8 {
        match self {
            FileBase::Bin => 0,
            FileBase::Resource => 1,
            FileBase::Share => 2,
        }
    }
}

/// Writes a streaming `NYXPK001` archive to `out_path` and returns TOC/hash metadata.
pub fn write_nyxpack(spec: &PackSpec, out_path: &Path) -> Result<WrittenPack> {
    if let Some(parent) = out_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp_body = out_path.with_extension(format!("body.{}.tmp", std::process::id()));
    let mut entries = Vec::with_capacity(spec.items.len());
    let mut data_offset: u64 = 0;

    {
        let body_file = File::create(&tmp_body)?;
        let mut body_writer = BufWriter::new(body_file);
        for item in &spec.items {
            validate_rel_path(&item.rel_path)?;
            let mut src = BufReader::new(
                File::open(&item.source_path)
                    .with_context(|| format!("opening '{}'", item.source_path.display()))?,
            );
            let mut hasher = blake3::Hasher::new();
            let mut encoder = ZlibEncoder::new(Vec::new(), Compression::new(6));
            let mut buf = [0u8; 64 * 1024];
            let mut raw_size = 0u64;
            loop {
                let n = src.read(&mut buf)?;
                if n == 0 {
                    break;
                }
                raw_size += n as u64;
                if raw_size > MAX_SINGLE_FILE_BYTES {
                    let _ = std::fs::remove_file(&tmp_body);
                    bail!("{} exceeds single-file cap", item.rel_path);
                }
                hasher.update(&buf[..n]);
                encoder.write_all(&buf[..n])?;
            }
            let compressed = encoder.finish()?;
            let comp_len = compressed.len() as u64;
            body_writer.write_all(&compressed)?;
            entries.push(BundleEntry {
                rel_path: item.rel_path.clone(),
                base: item.base,
                role: item.role,
                uncompressed_size: raw_size,
                compressed_size: comp_len,
                data_offset,
                blake3: format!("b3:{}", hasher.finalize().to_hex()),
                mode: item.mode,
            });
            data_offset += comp_len;
        }
        body_writer.flush()?;
    }

    let toc = BundleToc {
        format: 1,
        nyx_version: spec.nyx_version.clone(),
        target_triple: spec.target_triple.clone(),
        os: spec.os,
        entries,
    };
    let toc_bytes = serde_json::to_vec(&toc)?;
    let toc_len = u32::try_from(toc_bytes.len()).context("TOC exceeds 4 GiB")?;
    let toc_blake3 = *blake3::hash(&toc_bytes).as_bytes();

    let mut payload_hasher = blake3::Hasher::new();
    let mut out = BufWriter::new(File::create(out_path)?);
    let mut write_hashed = |slice: &[u8]| -> Result<()> {
        payload_hasher.update(slice);
        out.write_all(slice)?;
        Ok(())
    };
    write_hashed(PACK_MAGIC)?;
    write_hashed(&toc_len.to_le_bytes())?;
    write_hashed(&toc_bytes)?;

    let mut body_reader = BufReader::new(File::open(&tmp_body)?);
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = body_reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        write_hashed(&buf[..n])?;
    }
    drop(body_reader);
    let _ = std::fs::remove_file(&tmp_body);
    out.flush()?;

    let payload_len = 8 + 4 + (toc_bytes.len() as u64) + data_offset;
    let payload_blake3 = *payload_hasher.finalize().as_bytes();
    Ok(WrittenPack {
        toc_bytes,
        payload_len,
        toc_blake3,
        payload_blake3,
    })
}

/// Builds a self-extracting executable by cloning `stub_exe` (unsigned or pre-packed bare stub) and
/// appending `[NYXPK001 payload][0..7 pad bytes][96B NYXINST1 trailer]` aligned to an 8-byte boundary.
pub fn write_self_extracting(
    stub_exe: &Path,
    nyxpack_path: &Path,
    out_exe: &Path,
    toc_bytes: &[u8],
    toc_blake3: [u8; 32],
    payload_blake3: [u8; 32],
) -> Result<Trailer> {
    if let Some(parent) = out_exe.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp = out_exe.with_extension(format!("sfx.{}.tmp", std::process::id()));
    let mut out = BufWriter::new(File::create(&tmp)?);
    let stub_len = std::io::copy(&mut BufReader::new(File::open(stub_exe)?), &mut out)?;
    let payload_len = std::io::copy(&mut BufReader::new(File::open(nyxpack_path)?), &mut out)?;

    let after_payload = stub_len + payload_len;
    let pad = ((8 - ((after_payload + TRAILER_LEN as u64) % 8)) % 8) as usize;
    if pad > 0 {
        out.write_all(&vec![0u8; pad])?;
    }
    let trailer = Trailer {
        version: TRAILER_VERSION,
        codec: CODEC_ZLIB_RS,
        flags: 0,
        toc_len: u32::try_from(toc_bytes.len())?,
        payload_offset: stub_len,
        payload_len,
        toc_blake3,
        payload_blake3,
    };
    out.write_all(&trailer.to_bytes())?;
    out.flush()?;
    drop(out);
    if out_exe.exists() {
        let _ = std::fs::remove_file(out_exe);
    }
    std::fs::rename(&tmp, out_exe)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let _ = std::fs::set_permissions(out_exe, std::fs::Permissions::from_mode(0o755));
    }
    Ok(trailer)
}

/// Writes a standard `.zip` distribution archive containing the staged layout cleanly.
pub fn write_zip_bundle(spec: &PackSpec, out_zip: &Path) -> Result<()> {
    if let Some(parent) = out_zip.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let f = File::create(out_zip)?;
    let mut zw = zip::ZipWriter::new(f);
    for item in &spec.items {
        validate_rel_path(&item.rel_path)?;
        let archive_path = match item.base {
            FileBase::Bin => format!("bin/{}", item.rel_path),
            FileBase::Resource | FileBase::Share => item.rel_path.clone(),
        };
        let opts = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated)
            .unix_permissions(item.mode);
        zw.start_file(archive_path, opts)?;
        let mut src = BufReader::new(File::open(&item.source_path)?);
        std::io::copy(&mut src, &mut zw)?;
    }
    zw.finish()?;
    Ok(())
}

/// Stages all files in a `PackSpec` directly into a target layout directory (for unpacked dist or dry builds).
pub fn write_dir_tree(spec: &PackSpec, root_dir: &Path) -> Result<()> {
    std::fs::create_dir_all(root_dir)?;
    for item in &spec.items {
        validate_rel_path(&item.rel_path)?;
        let dst = match item.base {
            FileBase::Bin => root_dir.join("bin").join(&item.rel_path),
            FileBase::Resource => root_dir.join(&item.rel_path),
            FileBase::Share => root_dir.join("share/nyx").join(&item.rel_path),
        };
        if let Some(parent) = dst.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::copy(&item.source_path, &dst)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(item.mode));
        }
    }
    Ok(())
}

#[derive(Clone, Debug)]
pub struct ExtractedEntry {
    pub entry: BundleEntry,
    pub staged_path: PathBuf,
}

#[derive(Clone, Debug)]
pub struct OpenedPayload {
    pub toc: BundleToc,
    pub total_uncompressed: u64,
}

pub enum PayloadSource {
    Embedded { exe_path: PathBuf, trailer: Trailer },
    NyxPack { path: PathBuf },
    Zip { path: PathBuf },
    Dir { path: PathBuf },
    MacApp { app_path: PathBuf },
}

impl PayloadSource {
    /// Detects default source: `--bundle`, `--from-app`, `--from-dir`, embedded `NYXINST1`, or macOS sibling `../Nyx.app`.
    pub fn detect(
        current_exe: &Path,
        bundle: Option<&Path>,
        from_dir: Option<&Path>,
        from_app: Option<&Path>,
    ) -> Result<Option<Self>> {
        if let Some(p) = bundle {
            let ext = p
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            if ext == "zip" {
                return Ok(Some(Self::Zip {
                    path: p.to_path_buf(),
                }));
            }
            return Ok(Some(Self::NyxPack {
                path: p.to_path_buf(),
            }));
        }
        if let Some(app) = from_app {
            return Ok(Some(Self::MacApp {
                app_path: app.to_path_buf(),
            }));
        }
        if let Some(dir) = from_dir {
            return Ok(Some(Self::Dir {
                path: dir.to_path_buf(),
            }));
        }
        if let Ok(mut f) = File::open(current_exe) {
            if let Some((trailer, _)) = locate_embedded_trailer(&mut f)? {
                return Ok(Some(Self::Embedded {
                    exe_path: current_exe.to_path_buf(),
                    trailer,
                }));
            }
        }
        // Dual-path .dmg fallback: Set Up CLI & AI Hosts.app/Contents/MacOS/nyx-setup -> ../../../Nyx.app
        if let Some(app) = current_exe
            .canonicalize()
            .ok()
            .and_then(|p| p.ancestors().nth(4).map(|d| d.join("Nyx.app")))
            .filter(|a| a.join("Contents/MacOS/nyx").is_file())
        {
            return Ok(Some(Self::MacApp { app_path: app }));
        }
        Ok(None)
    }

    pub fn metadata_summary(&self, verifier: &dyn PayloadVerifier) -> Result<OpenedPayload> {
        match self {
            Self::Embedded { exe_path, trailer } => {
                let mut f = File::open(exe_path)?;
                read_nyxpack_toc(&mut f, trailer.payload_offset, Some(trailer), verifier)
            }
            Self::NyxPack { path } => {
                let mut f = File::open(path)?;
                read_nyxpack_toc(&mut f, 0, None, verifier)
            }
            Self::Zip { path } => inspect_zip(path),
            Self::Dir { path } => inspect_dir(path),
            Self::MacApp { app_path } => inspect_mac_app(app_path),
        }
    }

    /// Streams, decompresses, and BLAKE3-verifies all files into `stage_dir` inside `<base>/<rel_path>`.
    pub fn extract_to_stage(
        &self,
        stage_dir: &Path,
        target_os: TargetOs,
        uninstaller_stub: &Path,
        verifier: &dyn PayloadVerifier,
        on_entry: &mut dyn FnMut(usize, usize, &str),
    ) -> Result<BundleToc> {
        match self {
            Self::Embedded { exe_path, trailer } => {
                let mut f = BufReader::new(File::open(exe_path)?);
                extract_nyxpack(
                    &mut f,
                    trailer.payload_offset,
                    Some(trailer),
                    verifier,
                    stage_dir,
                    on_entry,
                )
            }
            Self::NyxPack { path } => {
                let mut f = BufReader::new(File::open(path)?);
                extract_nyxpack(&mut f, 0, None, verifier, stage_dir, on_entry)
            }
            Self::Zip { path } => extract_zip(path, target_os, stage_dir, on_entry),
            Self::Dir { path } => {
                copy_dir_source(path, target_os, uninstaller_stub, stage_dir, on_entry)
            }
            Self::MacApp { app_path } => {
                copy_mac_app_source(app_path, uninstaller_stub, stage_dir, on_entry)
            }
        }
    }
}

fn read_nyxpack_toc<R: Read + Seek>(
    r: &mut R,
    payload_offset: u64,
    trailer: Option<&Trailer>,
    verifier: &dyn PayloadVerifier,
) -> Result<OpenedPayload> {
    r.seek(SeekFrom::Start(payload_offset))?;
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic)?;
    if &magic != PACK_MAGIC {
        bail!("invalid payload magic (expected NYXPK001)");
    }
    let mut len_buf = [0u8; 4];
    r.read_exact(&mut len_buf)?;
    let toc_len = u32::from_le_bytes(len_buf);
    if toc_len > 16 * 1024 * 1024 {
        bail!("payload TOC exceeds 16 MiB limit");
    }
    let mut toc_bytes = vec![0u8; toc_len as usize];
    r.read_exact(&mut toc_bytes)?;
    if let Some(tr) = trailer {
        if tr.toc_len != toc_len {
            bail!("trailer TOC length mismatch");
        }
        verifier.verify_trailer(tr, &toc_bytes)?;
    }
    let toc: BundleToc = serde_json::from_slice(&toc_bytes)?;
    let mut total_uncompressed = 0u64;
    let mut seen = BTreeSet::new();
    for e in &toc.entries {
        validate_rel_path(&e.rel_path)?;
        if !seen.insert((e.base as u8, e.rel_path.clone())) {
            bail!("duplicate archive entry {}", e.rel_path);
        }
        if e.uncompressed_size > MAX_SINGLE_FILE_BYTES {
            bail!("archive entry {} exceeds size limit", e.rel_path);
        }
        total_uncompressed = total_uncompressed
            .checked_add(e.uncompressed_size)
            .context("total uncompressed overflow")?;
    }
    if total_uncompressed > MAX_TOTAL_PAYLOAD_BYTES {
        bail!("archive uncompressed size ({total_uncompressed}) exceeds 512 MiB cap");
    }
    Ok(OpenedPayload {
        toc,
        total_uncompressed,
    })
}

fn extract_nyxpack<R: Read + Seek>(
    r: &mut R,
    payload_offset: u64,
    trailer: Option<&Trailer>,
    verifier: &dyn PayloadVerifier,
    stage_dir: &Path,
    on_entry: &mut dyn FnMut(usize, usize, &str),
) -> Result<BundleToc> {
    let opened = read_nyxpack_toc(r, payload_offset, trailer, verifier)?;
    let data_start = payload_offset + 12;
    r.seek(SeekFrom::Start(payload_offset + 8))?;
    let mut lbuf = [0u8; 4];
    r.read_exact(&mut lbuf)?;
    let toc_len = u32::from_le_bytes(lbuf) as u64;
    let data_base = data_start + toc_len;

    if let Some(tr) = trailer {
        r.seek(SeekFrom::Start(payload_offset))?;
        let mut hasher = blake3::Hasher::new();
        let mut remaining = tr.payload_len;
        let mut buf = [0u8; 64 * 1024];
        while remaining > 0 {
            let chunk = (remaining as usize).min(buf.len());
            r.read_exact(&mut buf[..chunk])?;
            hasher.update(&buf[..chunk]);
            remaining -= chunk as u64;
        }
        verifier.verify_payload_stream(tr, hasher.finalize().as_bytes())?;
    }

    let total_entries = opened.toc.entries.len();
    let mut total_written = 0u64;
    for (idx, entry) in opened.toc.entries.iter().enumerate() {
        on_entry(idx + 1, total_entries, &entry.rel_path);
        let rel = validate_rel_path(&entry.rel_path)?;
        let base_dir = match entry.base {
            FileBase::Bin => stage_dir.join("bin"),
            FileBase::Resource => stage_dir.join("resource"),
            FileBase::Share => stage_dir.join("share"),
        };
        let dst = base_dir.join(rel);
        if let Some(p) = dst.parent() {
            std::fs::create_dir_all(p)?;
        }
        r.seek(SeekFrom::Start(data_base + entry.data_offset))?;
        let limited = r.by_ref().take(entry.compressed_size);
        let mut decoder = ZlibDecoder::new(limited);
        let mut out = BufWriter::new(File::create(&dst)?);
        let mut hasher = blake3::Hasher::new();
        let mut buf = [0u8; 64 * 1024];
        let mut raw_written = 0u64;
        loop {
            let n = decoder.read(&mut buf)?;
            if n == 0 {
                break;
            }
            raw_written += n as u64;
            total_written += n as u64;
            if raw_written > entry.uncompressed_size.max(MAX_SINGLE_FILE_BYTES)
                || total_written > MAX_TOTAL_PAYLOAD_BYTES
            {
                bail!("decompression limit exceeded for {}", entry.rel_path);
            }
            hasher.update(&buf[..n]);
            out.write_all(&buf[..n])?;
        }
        out.flush()?;
        if raw_written != entry.uncompressed_size {
            bail!(
                "size mismatch for {} (expected {}, got {raw_written})",
                entry.rel_path,
                entry.uncompressed_size
            );
        }
        let actual = format!("b3:{}", hasher.finalize().to_hex());
        if actual != entry.blake3 {
            bail!(
                "BLAKE3 checksum mismatch on extracted file {}",
                entry.rel_path
            );
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(entry.mode));
        }
    }
    Ok(opened.toc)
}

fn classify_zip_or_dir_rel(rel: &str, os: TargetOs) -> Result<(FileBase, String, FileRole, u32)> {
    let clean = rel.replace('\\', "/");
    let stripped = clean.strip_prefix("bin/").unwrap_or(&clean);
    let _ = validate_rel_path(stripped)?;
    let exe_suffix = if os == TargetOs::Windows { ".exe" } else { "" };
    if clean.starts_with("bin/") {
        let role = if stripped == format!("nyx-uninstall{exe_suffix}") {
            FileRole::Uninstaller
        } else {
            FileRole::Binary
        };
        return Ok((FileBase::Bin, stripped.to_string(), role, 0o755));
    }
    if matches!(
        clean.as_str(),
        "nyx"
            | "nyx.exe"
            | "nyx-mcp"
            | "nyx-mcp.exe"
            | "nyx-worker"
            | "nyx-worker.exe"
            | "nyx-script-worker"
            | "nyx-script-worker.exe"
            | "nyx-connect"
            | "nyx-connect.exe"
            | "nyx-uninstall"
            | "nyx-uninstall.exe"
    ) {
        let role = if clean.starts_with("nyx-uninstall") {
            FileRole::Uninstaller
        } else {
            FileRole::Binary
        };
        return Ok((FileBase::Bin, clean, role, 0o755));
    }
    let role = if matches!(
        clean.as_str(),
        "LICENSE" | "NOTICE" | "THIRD_PARTY_NOTICES.txt"
    ) {
        FileRole::License
    } else {
        FileRole::Asset
    };
    Ok((FileBase::Resource, clean, role, 0o644))
}

fn inspect_zip(path: &Path) -> Result<OpenedPayload> {
    let mut za = zip::ZipArchive::new(BufReader::new(File::open(path)?))?;
    let os = TargetOs::current();
    let mut entries = Vec::new();
    let mut total = 0u64;
    for i in 0..za.len() {
        let f = za.by_index(i)?;
        if f.is_dir() {
            continue;
        }
        let name = f.name()?.into_owned();
        if f.is_symlink() {
            bail!("zip archive contains unsupported symlink '{name}'");
        }
        if f.encrypted() {
            bail!("zip archive contains encrypted entry '{name}'");
        }
        let (base, rel_path, role, mode) = classify_zip_or_dir_rel(&name, os)?;
        if f.size() > MAX_SINGLE_FILE_BYTES {
            bail!("zip entry '{name}' exceeds size limit");
        }
        total = total.checked_add(f.size()).context("zip total overflow")?;
        if total > MAX_TOTAL_PAYLOAD_BYTES {
            bail!("zip uncompressed total exceeds 512 MiB limit");
        }
        entries.push(BundleEntry {
            rel_path,
            base,
            role,
            uncompressed_size: f.size(),
            compressed_size: f.compressed_size(),
            data_offset: 0,
            blake3: String::new(),
            mode,
        });
    }
    Ok(OpenedPayload {
        toc: BundleToc {
            format: 1,
            nyx_version: env!("CARGO_PKG_VERSION").to_string(),
            target_triple: target_triple_guess().to_string(),
            os,
            entries,
        },
        total_uncompressed: total,
    })
}

fn extract_zip(
    path: &Path,
    os: TargetOs,
    stage_dir: &Path,
    on_entry: &mut dyn FnMut(usize, usize, &str),
) -> Result<BundleToc> {
    let mut za = zip::ZipArchive::new(BufReader::new(File::open(path)?))?;
    let count = za.len();
    let mut entries = Vec::new();
    let mut total_written = 0u64;
    for i in 0..count {
        let mut f = za.by_index(i)?;
        if f.is_dir() {
            continue;
        }
        let name = f.name()?.into_owned();
        if f.is_symlink() {
            bail!("zip contains symlink entry '{name}'");
        }
        if f.encrypted() {
            bail!("zip contains encrypted entry '{name}'");
        }
        let (base, rel_path, role, mode) = classify_zip_or_dir_rel(&name, os)?;
        on_entry(i + 1, count, &rel_path);
        let sub = match base {
            FileBase::Bin => stage_dir.join("bin"),
            FileBase::Resource => stage_dir.join("resource"),
            FileBase::Share => stage_dir.join("share"),
        };
        let dst = sub.join(validate_rel_path(&rel_path)?);
        if let Some(p) = dst.parent() {
            std::fs::create_dir_all(p)?;
        }
        let mut hasher = blake3::Hasher::new();
        let mut out = BufWriter::new(File::create(&dst)?);
        let mut buf = [0u8; 64 * 1024];
        let mut raw = 0u64;
        loop {
            let n = f.read(&mut buf)?;
            if n == 0 {
                break;
            }
            raw += n as u64;
            total_written += n as u64;
            if raw > MAX_SINGLE_FILE_BYTES || total_written > MAX_TOTAL_PAYLOAD_BYTES {
                bail!("zip entry '{}' exceeds limit", rel_path);
            }
            hasher.update(&buf[..n]);
            out.write_all(&buf[..n])?;
        }
        out.flush()?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(mode));
        }
        entries.push(BundleEntry {
            rel_path,
            base,
            role,
            uncompressed_size: raw,
            compressed_size: f.compressed_size(),
            data_offset: 0,
            blake3: format!("b3:{}", hasher.finalize().to_hex()),
            mode,
        });
    }
    Ok(BundleToc {
        format: 1,
        nyx_version: env!("CARGO_PKG_VERSION").to_string(),
        target_triple: target_triple_guess().to_string(),
        os,
        entries,
    })
}

fn walk_dir_files(base: &Path, cur: &Path, out: &mut Vec<(String, PathBuf)>) -> Result<()> {
    let mut kids: Vec<PathBuf> = std::fs::read_dir(cur)?
        .filter_map(Result::ok)
        .map(|e| e.path())
        .collect();
    kids.sort();
    for p in kids {
        let meta = std::fs::symlink_metadata(&p)?;
        if meta.file_type().is_symlink() {
            continue;
        }
        if meta.is_dir() {
            walk_dir_files(base, &p, out)?;
        } else if meta.is_file() {
            let rel = p.strip_prefix(base)?.to_string_lossy().replace('\\', "/");
            if rel != "install-manifest.json" {
                out.push((rel, p));
            }
        }
    }
    Ok(())
}

fn inspect_dir(dir: &Path) -> Result<OpenedPayload> {
    if !dir.is_dir() {
        bail!("source directory does not exist: {}", dir.display());
    }
    let os = TargetOs::current();
    let mut raw_files = Vec::new();
    walk_dir_files(dir, dir, &mut raw_files)?;
    let mut entries = Vec::new();
    let mut total = 0u64;
    for (rel, path) in raw_files {
        let (base, rel_path, role, mode) = classify_zip_or_dir_rel(&rel, os)?;
        let size = std::fs::metadata(&path)?.len();
        total += size;
        entries.push(BundleEntry {
            rel_path,
            base,
            role,
            uncompressed_size: size,
            compressed_size: size,
            data_offset: 0,
            blake3: String::new(),
            mode,
        });
    }
    Ok(OpenedPayload {
        toc: BundleToc {
            format: 1,
            nyx_version: env!("CARGO_PKG_VERSION").to_string(),
            target_triple: target_triple_guess().to_string(),
            os,
            entries,
        },
        total_uncompressed: total,
    })
}

fn copy_dir_source(
    dir: &Path,
    os: TargetOs,
    uninstaller_stub: &Path,
    stage_dir: &Path,
    on_entry: &mut dyn FnMut(usize, usize, &str),
) -> Result<BundleToc> {
    if !dir.is_dir() {
        bail!("source directory does not exist: {}", dir.display());
    }
    let mut raw_files = Vec::new();
    walk_dir_files(dir, dir, &mut raw_files)?;
    let mut entries = Vec::new();
    let mut has_uninstaller = false;
    let total = raw_files.len();
    for (idx, (rel, src)) in raw_files.into_iter().enumerate() {
        let (base, rel_path, role, mode) = classify_zip_or_dir_rel(&rel, os)?;
        if role == FileRole::Uninstaller {
            has_uninstaller = true;
        }
        on_entry(idx + 1, total, &rel_path);
        let sub = match base {
            FileBase::Bin => stage_dir.join("bin"),
            FileBase::Resource => stage_dir.join("resource"),
            FileBase::Share => stage_dir.join("share"),
        };
        let dst = sub.join(validate_rel_path(&rel_path)?);
        if let Some(p) = dst.parent() {
            std::fs::create_dir_all(p)?;
        }
        let bytes = std::fs::read(&src)?;
        std::fs::write(&dst, &bytes)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(mode));
        }
        entries.push(BundleEntry {
            rel_path,
            base,
            role,
            uncompressed_size: bytes.len() as u64,
            compressed_size: bytes.len() as u64,
            data_offset: 0,
            blake3: format!("b3:{}", blake3::hash(&bytes).to_hex()),
            mode,
        });
    }
    if !has_uninstaller && uninstaller_stub.is_file() {
        let fname = exe_filename("nyx-uninstall", os);
        let dst = stage_dir.join("bin").join(&fname);
        if let Some(p) = dst.parent() {
            std::fs::create_dir_all(p)?;
        }
        let bytes = std::fs::read(uninstaller_stub)?;
        std::fs::write(&dst, &bytes)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(0o755));
        }
        entries.push(BundleEntry {
            rel_path: fname,
            base: FileBase::Bin,
            role: FileRole::Uninstaller,
            uncompressed_size: bytes.len() as u64,
            compressed_size: bytes.len() as u64,
            data_offset: 0,
            blake3: format!("b3:{}", blake3::hash(&bytes).to_hex()),
            mode: 0o755,
        });
    }
    Ok(BundleToc {
        format: 1,
        nyx_version: env!("CARGO_PKG_VERSION").to_string(),
        target_triple: target_triple_guess().to_string(),
        os,
        entries,
    })
}

fn inspect_mac_app(app_path: &Path) -> Result<OpenedPayload> {
    let contents = app_path.join("Contents");
    if !contents.is_dir() {
        bail!(
            "missing Contents/ inside macOS app bundle '{}'",
            app_path.display()
        );
    }
    let mut raw_files = Vec::new();
    walk_dir_files(&contents, &contents, &mut raw_files)?;
    let mut entries = Vec::new();
    let mut total = 0u64;
    for (rel, path) in raw_files {
        let size = std::fs::metadata(&path)?.len();
        total += size;
        let (base, rel_path, role, mode) = classify_mac_contents_rel(&rel)?;
        entries.push(BundleEntry {
            rel_path,
            base,
            role,
            uncompressed_size: size,
            compressed_size: size,
            data_offset: 0,
            blake3: String::new(),
            mode,
        });
    }
    Ok(OpenedPayload {
        toc: BundleToc {
            format: 1,
            nyx_version: env!("CARGO_PKG_VERSION").to_string(),
            target_triple: target_triple_guess().to_string(),
            os: TargetOs::Macos,
            entries,
        },
        total_uncompressed: total,
    })
}

fn classify_mac_contents_rel(rel: &str) -> Result<(FileBase, String, FileRole, u32)> {
    if let Some(rest) = rel.strip_prefix("MacOS/") {
        let role = if rest == "nyx-uninstall" {
            FileRole::Uninstaller
        } else {
            FileRole::Binary
        };
        return Ok((FileBase::Bin, rest.to_string(), role, 0o755));
    }
    if let Some(rest) = rel.strip_prefix("Resources/") {
        let role = if matches!(
            rest,
            "LICENSE" | "NOTICE" | "THIRD_PARTY_NOTICES.txt"
        ) {
            FileRole::License
        } else {
            FileRole::Asset
        };
        return Ok((FileBase::Resource, rest.to_string(), role, 0o644));
    }
    // Preserve Contents/Info.plist and Contents/_CodeSignature/* via FileBase::Share (mapped to Contents/).
    Ok((FileBase::Share, rel.to_string(), FileRole::Asset, 0o644))
}

fn copy_mac_app_source(
    app_path: &Path,
    uninstaller_stub: &Path,
    stage_dir: &Path,
    on_entry: &mut dyn FnMut(usize, usize, &str),
) -> Result<BundleToc> {
    let contents = app_path.join("Contents");
    if !contents.is_dir() {
        bail!(
            "missing Contents/ inside macOS app bundle '{}'",
            app_path.display()
        );
    }
    let mut raw_files = Vec::new();
    walk_dir_files(&contents, &contents, &mut raw_files)?;
    let mut entries = Vec::new();
    let mut has_uninstaller = false;
    let total = raw_files.len();
    for (idx, (rel, src)) in raw_files.into_iter().enumerate() {
        let (base, rel_path, role, mode) = classify_mac_contents_rel(&rel)?;
        if role == FileRole::Uninstaller {
            has_uninstaller = true;
        }
        on_entry(idx + 1, total, &rel_path);
        let sub = match base {
            FileBase::Bin => stage_dir.join("bin"),
            FileBase::Resource => stage_dir.join("resource"),
            FileBase::Share => stage_dir.join("share"),
        };
        let dst = sub.join(validate_rel_path(&rel_path)?);
        if let Some(p) = dst.parent() {
            std::fs::create_dir_all(p)?;
        }
        let bytes = std::fs::read(&src)?;
        std::fs::write(&dst, &bytes)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(mode));
        }
        entries.push(BundleEntry {
            rel_path,
            base,
            role,
            uncompressed_size: bytes.len() as u64,
            compressed_size: bytes.len() as u64,
            data_offset: 0,
            blake3: format!("b3:{}", blake3::hash(&bytes).to_hex()),
            mode,
        });
    }
    if !has_uninstaller && !contents.join("_CodeSignature").exists() && uninstaller_stub.is_file() {
        let dst = stage_dir.join("bin").join("nyx-uninstall");
        if let Some(p) = dst.parent() {
            std::fs::create_dir_all(p)?;
        }
        let bytes = std::fs::read(uninstaller_stub)?;
        std::fs::write(&dst, &bytes)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(0o755));
        }
        entries.push(BundleEntry {
            rel_path: "nyx-uninstall".to_string(),
            base: FileBase::Bin,
            role: FileRole::Uninstaller,
            uncompressed_size: bytes.len() as u64,
            compressed_size: bytes.len() as u64,
            data_offset: 0,
            blake3: format!("b3:{}", blake3::hash(&bytes).to_hex()),
            mode: 0o755,
        });
    }
    Ok(BundleToc {
        format: 1,
        nyx_version: env!("CARGO_PKG_VERSION").to_string(),
        target_triple: target_triple_guess().to_string(),
        os: TargetOs::Macos,
        entries,
    })
}

pub fn target_triple_guess() -> &'static str {
    if cfg!(all(target_os = "windows", target_arch = "x86_64")) {
        "x86_64-pc-windows-msvc"
    } else if cfg!(all(target_os = "windows", target_arch = "aarch64")) {
        "aarch64-pc-windows-msvc"
    } else if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
        "x86_64-unknown-linux-gnu"
    } else if cfg!(all(target_os = "linux", target_arch = "aarch64")) {
        "aarch64-unknown-linux-gnu"
    } else if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
        "aarch64-apple-darwin"
    } else if cfg!(all(target_os = "macos", target_arch = "x86_64")) {
        "x86_64-apple-darwin"
    } else {
        "unknown"
    }
}
