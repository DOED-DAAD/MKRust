use clap::{Parser, Subcommand};

mod fastq_quality;
mod kraken2;
mod mash;
mod mask;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None, arg_required_else_help=true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Mask(mask::MaskArgs),
    FastqQuality(fastq_quality::FastqQualityArgs),
    MashTopHit(mash::MashTopHitArgs),
    MashMetagenomic(mash::MashMetagenomicArgs),
    Kraken2TopHit(kraken2::Kraken2TopHitArgs),
    Kraken2Bin(kraken2::Kraken2BinArgs),
}

fn main() -> Result<(), std::io::Error> {
    let cli = Cli::parse();
    let _ = match cli.command {
        Commands::Mask(args) => mask::mask(args),
        Commands::FastqQuality(args) => fastq_quality::FastQQuality::get_fastq_stats(args),
        Commands::MashTopHit(args) => mash::get_mash_top_hit(args),
        Commands::MashMetagenomic(args) => mash::check_metagenomic(args),
        Commands::Kraken2TopHit(args) => kraken2::get_kraken2_top_hit(args),
        Commands::Kraken2Bin(args) => kraken2::bin_kraken2_fasta(args),
    };

    Ok(())
}
