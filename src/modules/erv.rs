//! Retroviral read detection module
//!
//! Screens reads for retroviral / EVE k-mer content and flags them
//! (informational; does not remove reads). The endogenous-vs-exogenous
//! classifier that previously consumed these flags was decommissioned
//! (see EVE_DISCRIMINATION_REDESIGN.md); the alignment-based EVE screen
//! (minimap2 + DIAMOND competitive EVE-vs-exogenous) is its replacement.

use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicU64, Ordering};

/// K-mer size for retroviral screening
const ERV_K: usize = 21;

/// Minimum k-mer containment fraction to flag a read as retroviral
const MIN_RETROVIRAL_FRACTION: f64 = 0.15;

/// Retroviral reference sequences compiled into the binary.
/// 44 Retroviridae genomes from RefSeq + HHV-6A/B for ciHHV-6 detection.
#[path = "erv_sequences.rs"]
pub(crate) mod erv_sequences;

/// Non-retroviral EVE reference sequences (Bornaviridae, Parvoviridae, Filoviridae).
#[path = "eve_sequences_nonretro.rs"]
pub(crate) mod eve_sequences;

/// EVE/ERV analysis module -- screens reads for viral k-mer content
///
/// Unlike other QcModule implementations, this module does NOT remove reads.
/// It flags reads with viral signal for downstream classification.
/// Covers Retroviridae (original) plus Bornaviridae, Parvoviridae, Filoviridae.
pub struct ErvScreener {
    /// Unified k-mer index (all viral families)
    retroviral_index: FxHashSet<u64>,
    /// Per-family k-mer sub-indices for family assignment
    family_indices: Vec<(&'static str, FxHashSet<u64>)>,
    /// Stats
    stats: AtomicStats,
    /// Reads with retroviral/viral signal
    retroviral_reads: AtomicU64,
}

impl Default for ErvScreener {
    fn default() -> Self {
        Self::new()
    }
}

/// Shared k-mer hashes between ERV retroviral index and herpesvirus genomes.
/// Computed from 30 herpesvirus genomes (HHV-1 through HHV-8 + common animal
/// herpesviruses). Removed from the ERV screening index to prevent false-positive
/// retroviral flagging on herpesvirus reads (DNA polymerase / RT domain homology).
/// 194 shared hashes.
const ERV_HERPES_SHARED: &[u64] = &[
    13962645681653164222, 4848821342062053254, 15254468161006421820,
    16720133836215698436, 11010584888746014358, 12825397808992046767,
    17385506404350697990, 4895285437823292849, 12924188559524224447,
    65195903065308832, 10081571324861535705, 1049834778505101715,
    18363656723753895607, 16727997448078670466, 15851609270537664653,
    3071281897178906356, 2464970366815453046, 9940973495366111066,
    13319355644564361605, 17932056254999167885, 10565343728905348550,
    6046725081124947622, 6046741573799370787, 4000368309555412237,
    6857969769562877030, 18139156394498853206, 12299362826762190545,
    3058439901316153929, 10870943296901293006, 6977631160319285446,
    12303828737447554347, 15071834637225634644, 17527370655297203544,
    10970495850922510826, 12883888612314180798, 2193259526192046883,
    10081590016559215292, 3482159967714402531, 8510924157394611038,
    10086475316878740575, 17549325623098315601, 17107210390128727378,
    9914607167265701131, 3069102965084462732, 18201028470742597975,
    15081171782337311446, 12174930535191439862, 9739535805503296490,
    15758046215814736675, 11804167051853324063, 17229401272156212687,
    923452153065247026, 11043983704436214096, 1941470522670562861,
    7162062784630449948, 17937892462720576523, 18157815932969661759,
    18430728532443968368, 6977633359342541868, 14945647260567753618,
    6055183624079071820, 17259671278210084469, 1844515583590083098,
    5102305008234515306, 9597352542161744647, 10066252749985580471,
    6241116355873728466, 1819474666958779932, 13993129839198383702,
    3877902065016142710, 4227142601420069988, 10282317800072864950,
    3692604543722650914, 14683054609334607875, 5792920969797081691,
    5390957747897509612, 14764903304866444388, 8810438012324563156,
    6046729479171460466, 13030267304690332738, 2344286985023316642,
    2073487729902070973, 12299367224808703389, 6857990660283813039,
    17138914518006727951, 652487057309969677, 1336345004851995159,
    12194612728323611720, 7016202814331052675, 8161033317928423250,
    103972683079733980, 342153770148907437, 11784535694598692678,
    7922478350982330354, 7713930007535559067, 3451449663956728333,
    13266146183837845441, 12681940554868902383, 1523061571148126771,
    5879471276037682154, 11043988102482726940, 12479372993012529594,
    18044759732823680454, 526072545103210290, 11467371894108803459,
    3070873178805503992, 1899870825726972586, 3487826069628704090,
    5877429482944472777, 16363539529624384292, 681610932426196060,
    8382050345048132193, 9551344220846356531, 7000759313081216694,
    2922837332283829660, 18229992827549127631, 14943986846912971848,
    6966892652686386484, 12877684179448114099, 3747758593027569358,
    7901272927439122794, 16681927381889026614, 11060816371201573613,
    1026150274276257978, 11043990301505983362, 9599270090440966181,
    6586363147478644533, 9782844556924593718, 15581838239960184278,
    14953139160715486848, 9417444961663044841, 15276448497961635746,
    3070875377828760414, 13997860522397084794, 2236385220494696304,
    17087802327149262358, 5461880994830361015, 11475813944388504492,
    8894470780115434753, 3760055531075021257, 1664226575132881862,
    11912838100681596186, 8651969594106463900, 3749442038141957229,
    2731433690355295990, 4571699415554723727, 16630854718614439100,
    2616245650830005782, 7713917912907648746, 16664234982446853356,
    6386338672418595996, 15342774367820847184, 11144496540758692422,
    16681929580912283036, 13847057871564618086, 10767201891715687419,
    5375022635640163092, 6059162756660810529, 9687434104115470107,
    7085616550517903010, 9782846755947850140, 17350353325426591902,
    10029570010857037303, 12175102611877824258, 4263478609198614115,
    11010566197048334771, 9571956392821920742, 6829646631662399907,
    11535041740783142495, 16630856917637695522, 6933276945136837883,
    7851875135511764119, 3733178062140676942, 16565651194698782422,
    6386324378767429253, 12695861413934518864, 6761700589962934710,
    9068681361804278023, 2276580127378307809, 11584600435737316561,
    3890357974768973590, 16851302816060397641, 3947883906524771171,
    8099888696268154806, 14415517647924486404, 11629168612341181107,
    17761334120765046355, 3243714519403542382, 3071279698155649934,
    15814034966919130738, 2803395663069079934, 5086975889742534608,
    17855744689538449507, 8209570544467295562,
];

impl ErvScreener {
    /// Build the unified viral k-mer index from all embedded sequences.
    /// Includes Retroviridae (original) plus Bornaviridae, Parvoviridae, Filoviridae.
    pub fn new() -> Self {
        let mut index = FxHashSet::default();
        let mut family_indices: Vec<(&str, FxHashSet<u64>)> = Vec::new();

        // Retroviridae (existing)
        let mut retro_idx = FxHashSet::default();
        for seq in erv_sequences::RETROVIRAL_SEQUENCES {
            add_sequence_kmers(seq, ERV_K, &mut index);
            add_sequence_kmers(seq, ERV_K, &mut retro_idx);
        }
        family_indices.push(("Retroviridae", retro_idx));

        // Remove k-mers shared with herpesvirus genomes to prevent
        // false-positive retroviral flagging on herpesvirus reads
        let before = index.len();
        for &shared_hash in ERV_HERPES_SHARED {
            index.remove(&shared_hash);
        }
        let herpes_removed = before - index.len();

        // Non-retroviral EVE families
        for &(name, seq) in eve_sequences::ALL_EVE_SEQUENCES {
            let family = name.split(':').next().unwrap_or("unknown");
            add_sequence_kmers(seq, ERV_K, &mut index);

            // Add to family-specific index
            if let Some((_, ref mut fam_idx)) = family_indices.iter_mut().find(|(f, _)| *f == family) {
                add_sequence_kmers(seq, ERV_K, fam_idx);
            } else {
                let mut fam_idx = FxHashSet::default();
                add_sequence_kmers(seq, ERV_K, &mut fam_idx);
                family_indices.push((
                    // Leak the family string since we need 'static lifetime
                    // Safe because these are compile-time constants from the embedded sequences
                    Box::leak(family.to_string().into_boxed_str()),
                    fam_idx,
                ));
            }
        }

        let eve_families: Vec<String> = family_indices
            .iter()
            .map(|(name, idx)| format!("{} ({})", name, idx.len()))
            .collect();

        log::info!(
            "EVE/ERV screening index: {} unique k-mers ({} herpes-shared removed). Families: {}",
            index.len(),
            herpes_removed,
            eve_families.join(", "),
        );

        Self {
            retroviral_index: index,
            family_indices,
            stats: AtomicStats::new(),
            retroviral_reads: AtomicU64::new(0),
        }
    }

    /// Determine which viral family a read most likely belongs to
    /// based on k-mer hits against family-specific sub-indices.
    /// Returns the family name with the highest containment.
    pub fn dominant_family(&self, sequence: &[u8]) -> Option<&'static str> {
        if sequence.len() < ERV_K || self.family_indices.is_empty() {
            return None;
        }

        let mut best_family = None;
        let mut best_hits = 0u32;

        for &(family_name, ref fam_idx) in &self.family_indices {
            let mut hits = 0u32;
            for kmer in sequence.windows(ERV_K) {
                if fam_idx.contains(&hash_kmer(kmer)) {
                    hits += 1;
                }
            }
            if hits > best_hits {
                best_hits = hits;
                best_family = Some(family_name);
            }
        }

        best_family
    }

    /// Compute retroviral k-mer containment of a read
    fn containment(&self, sequence: &[u8]) -> f64 {
        if sequence.len() < ERV_K || self.retroviral_index.is_empty() {
            return 0.0;
        }
        let total = sequence.len() - ERV_K + 1;
        let mut hits = 0u32;
        for kmer in sequence.windows(ERV_K) {
            if self.retroviral_index.contains(&hash_kmer(kmer)) {
                hits += 1;
            }
        }
        hits as f64 / total as f64
    }
}

impl QcModule for ErvScreener {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.len() < ERV_K {
            return;
        }

        let containment = self.containment(seq);
        if containment >= MIN_RETROVIRAL_FRACTION {
            // Record retroviral signal in metrics but do NOT change disposition.
            // The read stays in clean output -- ERV analysis is informational.
            // Using metrics rather than flag() to avoid routing to ambiguous output.
            record.metrics.retroviral_signal = true;
            self.retroviral_reads.fetch_add(1, Ordering::Relaxed);
        }
    }

    fn report(&self) -> ModuleReport {
        let retroviral = self.retroviral_reads.load(Ordering::Relaxed);
        let families: serde_json::Value = self
            .family_indices
            .iter()
            .map(|(name, idx)| (name.to_string(), serde_json::json!(idx.len())))
            .collect::<serde_json::Map<String, serde_json::Value>>()
            .into();

        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "retroviral_reads_flagged": retroviral,
                "index_size": self.retroviral_index.len(),
                "family_index_sizes": families,
            }),
        )
    }

    fn name(&self) -> &str {
        "erv"
    }
}

/// FNV-1a hash (uppercase, same as contaminant screener)
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Hash reverse complement
#[inline]
fn hash_kmer_rc(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        };
        hash ^= comp as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Add all k-mers from a sequence to a hash set (both orientations)
fn add_sequence_kmers(seq: &[u8], k: usize, set: &mut FxHashSet<u64>) {
    if seq.len() < k {
        return;
    }
    for kmer in seq.windows(k) {
        if kmer.iter().any(|&b| b == b'N' || b == b'n') {
            continue;
        }
        set.insert(hash_kmer(kmer));
        set.insert(hash_kmer_rc(kmer));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use biometal::FastqRecord;

    fn make_record(seq: &[u8]) -> AnnotatedRecord {
        let qual = vec![b'I'; seq.len()];
        AnnotatedRecord::new(FastqRecord::new("test".into(), seq.to_vec(), qual))
    }

    #[test]
    fn test_erv_index_construction() {
        let screener = ErvScreener::new();
        assert!(
            screener.retroviral_index.len() > 10_000,
            "Should have substantial retroviral k-mer index, got {}",
            screener.retroviral_index.len()
        );
    }

    #[test]
    fn test_retroviral_read_flagged() {
        let screener = ErvScreener::new();
        // Take a fragment from HIV-1 (first retroviral sequence)
        let hiv_fragment = &erv_sequences::RETROVIRAL_SEQUENCES[0][100..250];
        let mut record = make_record(hiv_fragment);
        screener.process(&mut record);
        // Should have retroviral_signal set (not flagged or failed)
        assert!(
            record.metrics.retroviral_signal,
            "Retroviral read should have retroviral_signal metric set"
        );
        assert!(
            !record.is_failed(),
            "Retroviral read should NOT be failed (ERV module is informational)"
        );
        assert!(
            !record.is_flagged(),
            "Retroviral read should NOT be flagged (stays in clean output)"
        );
    }

    #[test]
    fn test_random_read_not_flagged() {
        let screener = ErvScreener::new();
        let random: Vec<u8> = (0..150)
            .map(|i| [b'A', b'C', b'G', b'T'][(i * 7 + 3) % 4])
            .collect();
        let mut record = make_record(&random);
        screener.process(&mut record);
        assert!(!record.is_flagged(), "Random read should not be flagged");
    }
}
