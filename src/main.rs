use clap::clap_derive;
use clap::Parser;

use indicatif::MultiProgress;
use indicatif::ProgressBar;
use std::collections::HashMap;
use std::io::BufRead;
use std::io::Write;
use std::path::PathBuf;

use md5;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::thread;

#[derive(clap_derive::Parser, Debug)]
#[clap(
    name = "smart_gnomad_downloader",
    version = "1.0",
    author = "SteampunkIslande"
)]
#[clap(about = "Download bed-restricted gnomad genome vcf")]
struct Cli {
    /// Path to the BED file to restrict the VCF to
    #[clap(short = 'b', long = "bed")]
    bed: PathBuf,

    /// Path to the BED file to URL CSV list
    #[clap(short = 'u', long = "url-list")]
    urls: PathBuf,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct URLDownloadRecord {
    chromosome: String,
    md5sum: String,
    url: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct BEDRecord {
    chromosome: String,
    start: u32,
    end: u32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct VCFRecord {
    chromosome: String,
    pos: u32,
}

struct Md5ConsumerWriter {
    md5_context: md5::Context,
    progress_bar: ProgressBar,
}

impl Write for Md5ConsumerWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.md5_context.consume(buf);
        let chunk_length = buf.len();
        self.progress_bar.inc(chunk_length as u64);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.md5_context.flush()
    }
}

impl Md5ConsumerWriter {
    pub fn new(progress_bar: ProgressBar) -> Self {
        Md5ConsumerWriter {
            md5_context: md5::Context::new(),
            progress_bar,
        }
    }

    pub fn digest(self) -> String {
        self.progress_bar.finish();
        format!("{:x}", self.md5_context.compute())
    }
}

struct SortedIntervalIntersect<I, T>
where
    I: Iterator<Item = (T, T)>,
{
    intervals: I,
    current_interval: Option<(T, T)>,
}

impl<I, T> SortedIntervalIntersect<I, T>
where
    I: Iterator<Item = (T, T)>,
    T: PartialOrd,
{
    pub fn new(intervals: I) -> Self {
        Self {
            intervals,
            current_interval: None,
        }
    }

    pub fn in_interval(&mut self, value: T) -> Option<bool> {
        if self.current_interval.is_none() {
            self.current_interval = self.intervals.next();
        }
        loop {
            if let Some(ref current_interval) = self.current_interval {
                if value < current_interval.0 {
                    return Some(false);
                }
                if value < current_interval.1 {
                    return Some(true);
                }
                self.current_interval = self.intervals.next();
            } else {
                return None;
            }
        }
    }
}

fn get_blocking_reader_from_url(url: &str) -> Option<reqwest::blocking::Response> {
    let client = reqwest::blocking::Client::builder()
        .timeout(None)
        .build()
        .expect("Cannot build client");
    client.get(url).send().ok()
}

fn smart_save_vcf_from_url<I>(
    url: &str,
    expected_md5: &str,
    regions: I,
    output_file_name: &str,
    progress_bar: ProgressBar,
) where
    I: Iterator<Item = (u32, u32)>,
{
    eprintln!("Downloading file from {}", url);

    let raw_reader =
        get_blocking_reader_from_url(url).expect(&format!("Cannot download file at {}", url));

    let pb = progress_bar;
    if let Some(len) = raw_reader.content_length() {
        pb.set_length(len);
    }

    let mut md5_writer = Md5ConsumerWriter::new(pb);

    let actual_reader = tee::TeeReader::new(raw_reader, &mut md5_writer);
    let bg_reader = noodles::bgzf::Reader::new(actual_reader);

    let vcf_file =
        std::fs::File::create(&format!("{}", output_file_name)).expect("Cannot create output file");

    let mut vcf_file_writer = noodles::bgzf::io::Writer::new(vcf_file);

    let mut intersection_check = SortedIntervalIntersect::new(regions);

    for line in bg_reader
        .lines()
        .map(Result::ok)
        .map(|o| o.expect("Cannot read line!"))
    {
        if line.starts_with("#") {
            writeln!(vcf_file_writer, "{}", line)
                .expect("Cannot write to output! Received EOF or whatever...");
        } else {
            let pos: u32 = line
                .split("\t")
                .nth(1)
                .expect(&format!("Invalid vcf line: {}", line))
                .parse()
                .expect(&format!(
                    "{}: second field is not a number, or is larger than 4 billion!",
                    line
                ));
            if let Some(intersects) = intersection_check.in_interval(pos) {
                if intersects {
                    writeln!(vcf_file_writer, "{}", line)
                        .expect("Cannot write to output! Received EOF or whatever...");
                }
            }
        }
    }

    let downloaded_md5 = md5_writer.digest();
    if &downloaded_md5 == expected_md5 {
        eprintln!(
            "MD5 checksums match! The whole VCF was downloaded properly! - {}",
            output_file_name
        );
    } else {
        eprintln!(
            "MD5 checksums do not match! Part of the VCf was corrupted during download! - {}",
            output_file_name
        );
    }
}

fn main() {
    let args = Cli::parse();

    let bed_path = args.bed;
    let urls_path = args.urls;

    let mut bed_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(File::open(bed_path.as_path()).unwrap());

    let mut urls_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b',')
        .flexible(true)
        .from_reader(File::open(urls_path.as_path()).unwrap());

    let urls: HashMap<String, (String, String)> = urls_reader
        .deserialize::<URLDownloadRecord>()
        .into_iter()
        .map(|result| result.expect("Failed to deserialize URL record"))
        .fold(HashMap::new(), |mut acc, record| {
            acc.entry(record.chromosome)
                .or_insert((record.md5sum, record.url));
            acc
        });

    // Assumes the bed is not sorted. Could be optimized if the user could guarantee the BED is sorted...
    let regions_per_chr: HashMap<String, Vec<(u32, u32)>> = {
        let mut regions = bed_reader
            .deserialize::<BEDRecord>()
            .into_iter()
            .map(|result| result.expect("Failed to deserialize BED record"))
            .fold(HashMap::new(), |mut acc, record| {
                acc.entry(record.chromosome)
                    .or_insert_with(Vec::new)
                    .push((record.start, record.end));
                acc
            });

        for (_, regions) in regions.iter_mut() {
            regions.sort_by_key(|r| r.0);
        }
        regions
    };

    let mut thread_handles: Vec<thread::JoinHandle<_>> = Vec::new();

    let multi_progress: MultiProgress = MultiProgress::new();

    for (chrom_name, regions) in regions_per_chr.into_iter() {
        if let Some((expected_md5, url)) = urls.get(&chrom_name) {
            let url = url.clone();
            let expected_md5 = expected_md5.clone();
            let chrom_name = chrom_name.clone();
            let regions = regions.clone();

            let pb = multi_progress.add(ProgressBar::no_length());

            thread_handles.push(thread::spawn(move || {
                smart_save_vcf_from_url(
                    &url,
                    &expected_md5,
                    regions.into_iter(),
                    &format!("{}.vcf.gz", &chrom_name),
                    pb,
                )
            }));
        }
    }
    for handle in thread_handles {
        if let Err(e) = handle.join() {
            eprintln!("Thread panicked with message {:?}", e);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_intersection_check() {
        let regions = [(100u32, 200u32), (400u32, 1000u32)];
        let mut intersection_check = SortedIntervalIntersect::new(regions.into_iter());

        assert_eq!(intersection_check.in_interval(10), Some(false));
        assert_eq!(intersection_check.in_interval(150), Some(true));
        assert_eq!(intersection_check.in_interval(151), Some(true));
        assert_eq!(intersection_check.in_interval(220), Some(false));
        assert_eq!(intersection_check.in_interval(450), Some(true));
        assert_eq!(intersection_check.in_interval(460), Some(true));
        assert_eq!(intersection_check.in_interval(1100), None);
        assert_eq!(intersection_check.in_interval(2000), None);
        assert_eq!(intersection_check.in_interval(3000), None);
    }
}
