use clap::{Args, ValueEnum};
use flate2::Compression;
use flate2::write::GzEncoder;
use regex::Regex;
use rust_htslib::faidx::{Reader, build};

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

macro_rules! define_enum_with_strings {
    (
        $(#[$meta:meta])*
        $vis:vis enum $name:ident {
            $($variant:ident),* $(,)?
        }
    ) => {
        $(#[$meta])*
        $vis enum $name {
            $($variant),*
        }

        impl $name {


            //fn to_string(&self) -> &'static str {
            //    match self {
            //        $( $name::$variant => stringify!($variant),)*
            //    }
            //}

            fn from_string(s: &str) -> TaxonomicLevel {
                match s {
                    $( stringify!($variant) => $name::$variant,)*
                    _ => panic!("Invalid Taxonomic Level {}", s)
                }
            }

        }
    };
}

define_enum_with_strings!(
    #[derive(Copy, Clone, PartialEq, PartialOrd, Ord, Eq, ValueEnum, Debug, Hash)]
    enum TaxonomicLevel {
        S4,
        S3,
        S2,
        S1,
        S,
        G3,
        G2,
        G1,
        G,
        F3,
        F2,
        F1,
        F,
        O4,
        O3,
        O2,
        O1,
        O,
        C6,
        C5,
        C4,
        C3,
        C2,
        C1,
        C,
        P9,
        P8,
        P7,
        P6,
        P5,
        P4,
        P3,
        P2,
        P1,
        P,
        K,
        D5,
        D4,
        D3,
        D2,
        D1,
        D,
        R3,
        R2,
        R1,
        R,
        U,
    }
);

#[allow(dead_code)]
struct Kraken2ReportRow {
    percentage: f64,
    fragments_recovered: u64,
    fragments_assigned_taxon: u64,
    rank_code: TaxonomicLevel,
    taxon_id: u64,
    taxon_name: String,
}

impl Kraken2ReportRow {
    fn from_line(line: &str) -> Self {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        Self {
            percentage: fields[0]
                .parse::<f64>()
                .expect("Failed to parse percentage"),
            fragments_recovered: fields[1]
                .parse::<u64>()
                .expect("Failed to parse fragments recovered"),
            fragments_assigned_taxon: fields[2]
                .parse::<u64>()
                .expect("Failed to parse fragments assigned taxon"),
            rank_code: TaxonomicLevel::from_string(fields[3]),
            taxon_id: fields[4].parse::<u64>().expect("Failed to parse taxon ID"),
            taxon_name: fields[5].trim().to_string(),
        }
    }
}

#[derive(Args, Debug)]
pub struct Kraken2TopHitArgs {
    #[arg(short, long, help = "Input Kraken2 report")]
    report: PathBuf,

    #[arg(short, long, value_enum, default_value_t=TaxonomicLevel::S, help = "Taxonomic level to classify at.") ]
    taxa_level: TaxonomicLevel,
}

#[derive(Args, Debug)]
pub struct Kraken2BinArgs {
    #[arg(short, long, help = "Input Kraken2 report.")]
    report: PathBuf,

    #[arg(short, long, help = "Input Kraken2 output.")]
    kraken_output: PathBuf,

    #[arg(short, long, help = "Input fasta")]
    contigs: PathBuf,

    #[arg(short, long, help = "Output directory for binned fastas.")]
    prefix: String,

    #[arg(short, long, value_enum, default_value_t=TaxonomicLevel::S, help = "Taxonomic level to classify at.") ]
    taxa_level: TaxonomicLevel,

    #[arg(
        long,
        default_value_t = false,
        help = "Keep contigs classified at a higher level."
    )]
    keep_parents: bool,
}

fn get_taxa_id_sets(
    args: &Kraken2BinArgs,
) -> Result<HashMap<String, HashSet<u64>>, std::io::Error> {
    let report_path = args.report.clone();
    let file = File::open(report_path)?;
    let mut buf = String::new();
    let mut reader = BufReader::new(file);
    let mut taxa_stack: Vec<TaxonomicLevel> = Vec::new();
    let mut taxid_stack: Vec<u64> = Vec::new();
    let mut taxid_sets: HashMap<String, HashSet<u64>> = HashMap::new();
    let mut current_value: Option<&mut HashSet<u64>> = None;

    /*
     * Algo notes:
     * Make a Set containing our desired binning level and add taxids of all
     * children at that level and parent taxid if keep parents is selected.
     *
     * Then as we read the kraken2 output file we can just test if the taxid is in that
     * set then write it out.
     * */

    while reader.read_line(&mut buf).unwrap() > 0 {
        if buf.starts_with('#') {
            buf.clear();
            continue; // Skip comment lines
        }

        let kraken2_report_row = Kraken2ReportRow::from_line(&buf);

        if kraken2_report_row.rank_code == TaxonomicLevel::U {
            buf.clear();
            continue; // Skip unclassified rows
        }

        if let Some(last_rank) = taxa_stack.last()
            && *last_rank > kraken2_report_row.rank_code
        {
            taxa_stack.push(kraken2_report_row.rank_code);
            taxid_stack.push(kraken2_report_row.taxon_id);
        } else {
            while let Some(last_rank) = taxa_stack.last() {
                if *last_rank <= kraken2_report_row.rank_code {
                    taxa_stack.pop();
                    taxid_stack.pop();
                } else {
                    break;
                }
            }
            taxa_stack.push(kraken2_report_row.rank_code);
            taxid_stack.push(kraken2_report_row.taxon_id);
        }

        // TODO: Update docs to mention if the selected value is does not exist
        // for a value e.g. S2 for supbtypes. Than anything that does not have
        // a S2 lelvel will not show up.
        if args.taxa_level == kraken2_report_row.rank_code {
            // Now start filling set structure with tax ids
            taxid_sets
                .entry(kraken2_report_row.taxon_name.clone())
                .or_default()
                .insert(kraken2_report_row.taxon_id);

            if args.keep_parents {
                // Add the parent values to the set if needed
                taxid_stack.iter().for_each(|taxid| {
                    taxid_sets
                        .get_mut(&kraken2_report_row.taxon_name)
                        .unwrap()
                        .insert(*taxid);
                });
            }
            current_value = Some(taxid_sets.get_mut(&kraken2_report_row.taxon_name).unwrap());
        }

        if let Some(ref mut curr_set) = current_value
            && kraken2_report_row.rank_code < args.taxa_level
        {
            // Add the parent values to the set if needed
            curr_set.insert(kraken2_report_row.taxon_id);
        }

        buf.clear();
    }

    Ok(taxid_sets)
}

fn get_kraken2_contig_ids(args: &Kraken2BinArgs) -> Vec<(String, u64)> {
    let mut contig_tax_id: Vec<(String, u64)> = Vec::new();

    let file = File::open(args.kraken_output.clone());
    let mut buf = String::new();
    let mut reader = BufReader::new(file.unwrap());

    while reader.read_line(&mut buf).unwrap() > 0 {
        if buf.starts_with('#') {
            buf.clear();
            continue; // Skip comment lines
        }
        let fields: Vec<&str> = buf.trim().split('\t').collect();
        // Dont write out unclassified as that is in its own file from kraken2
        let classified = fields[0] == "C";
        if !classified {
            buf.clear();
            continue;
        }
        let contig_id = fields[1].to_string();
        let tax_id = fields[2].trim().parse::<u64>().unwrap();
        contig_tax_id.push((contig_id, tax_id));
        buf.clear();
    }
    contig_tax_id
}

fn get_contigs(
    args: &Kraken2BinArgs,
    ordered_contigs: Vec<(String, String)>,
) -> Result<(), std::io::Error> {
    let contigs_path = args.contigs.clone();
    let faidx_path = contigs_path.with_added_extension("fai");
    if !faidx_path.exists() {
        build(&contigs_path).expect("Failed to build fasta index.");
    }
    let mut output_path = std::env::current_dir()?;
    // Supports both compressed and uncompressed
    let fasta_reader = Reader::from_path(contigs_path).expect("Failed to open fasta file.");

    let mut first_output = ordered_contigs
        .first()
        .map(|(_, taxa_name)| taxa_name.clone())
        .expect("No classified contigs classified at selected output level.");
    let output_closure = |output_name| format!("{}_{}.binned.fasta.gz", args.prefix, output_name);
    output_path.push(output_closure(first_output.clone()));
    let mut file = File::create(output_path)?;
    let mut buf_writer = BufWriter::new(file);
    let mut encoder = GzEncoder::new(buf_writer, Compression::default());

    for (contig_name, taxa_name) in ordered_contigs {
        if taxa_name != first_output {
            buf_writer = encoder.finish()?;
            buf_writer.flush()?;
            first_output = taxa_name.clone();
            output_path = std::env::current_dir()?;
            output_path.push(output_closure(first_output.clone()));
            file = File::create(output_path)?;
            buf_writer = BufWriter::new(file);
            encoder = GzEncoder::new(buf_writer, Compression::default());
        }
        if let Ok(sequence) = fasta_reader.fetch_seq(
            &contig_name,
            0,
            fasta_reader.fetch_seq_len(&contig_name) as usize,
        ) {
            encoder.write_all(b">")?;
            encoder.write_all(contig_name.as_bytes())?;
            encoder.write_all(b"\n")?;
            encoder.write_all(&sequence)?;
            encoder.write_all(b"\n")?;
        }
    }
    buf_writer = encoder.finish()?;
    buf_writer.flush()?;
    Ok(())
}

pub fn bin_kraken2_fasta(args: Kraken2BinArgs) -> Result<(), std::io::Error> {
    let taxa_id_sets = get_taxa_id_sets(&args)?;
    let contig_tax_ids = get_kraken2_contig_ids(&args);
    // Need to process to make a look up for the vectors first
    //
    // Can make new hashmap where each contig string is the key, and we can look
    // up the taxid in the taxa_id_sets, and store the organism name along side it in a tuple.
    //
    // Need to order writes based on the fasta first, perhaps using faidx

    let re = Regex::new(r"[^A-Za-z0-9\-_]").unwrap();

    let mut matched_ids: Vec<(String, String)> = Vec::new();
    for (name, taxid) in &contig_tax_ids {
        for (taxa_name, taxa_set) in &taxa_id_sets {
            if taxa_set.contains(taxid) {
                // Iterate over full set if as parent contigs can appear more than once
                let replaced_text = re.replace_all(taxa_name, "_");
                matched_ids.push((name.clone(), replaced_text.to_string()));
            }
        }
    }

    matched_ids.sort_by(|a, b| a.1.cmp(&b.1));

    get_contigs(&args, matched_ids)?;

    Ok(())
}

pub fn get_kraken2_top_hit(args: Kraken2TopHitArgs) -> Result<(), std::io::Error> {
    let report_path = args.report;
    let taxa_level = args.taxa_level;

    let file = File::open(report_path)?;
    let mut reader = BufReader::new(file);

    let mut buf = String::new();

    let mut kraken2_top_hits: HashMap<TaxonomicLevel, Kraken2ReportRow> = HashMap::new();

    while reader.read_line(&mut buf).unwrap() > 0 {
        if buf.starts_with('#') {
            buf.clear();
            continue; // Skip comment lines
        }
        let kraken2_report_row = Kraken2ReportRow::from_line(&buf);

        match kraken2_top_hits.entry(kraken2_report_row.rank_code) {
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(kraken2_report_row);
            }
            // TODO: Need to add additional logic if the top-hit has more than one value
            std::collections::hash_map::Entry::Occupied(mut e) => {
                if kraken2_report_row.percentage > e.get().percentage {
                    e.insert(kraken2_report_row);
                }
            }
        }
        buf.clear();
    }

    let mut output_value = "No Species Identified".to_string();
    if let Some(top_hit) = kraken2_top_hits.get(&taxa_level) {
        output_value = top_hit.taxon_name.clone().replace('"', "");
    }
    print!("{}", output_value);

    Ok(())
}
