# Why smart_gnomad_downloader?

This simple rust project, with no dependencies, allows one to download big VCF files such as gnomad, while only saving to disk what intersects a given BED file.

I know you can simply `wget` or `curl` the file and pipe it to `bedtools`, but this can be tedious, especially if you need to ensure the checksum matches.

# How to install

Assuming you already have Rust installed:

```bash
git clone https://github.com/SteampunkIslande/smart_gnomad_downloader
cargo install --path .
```

Otherwise, open an issue and I will provide you with the binaries.

# How to use

```bash
smart_gnomad_downloader --bed test.bed --url-list gnomad_urls_download.csv
```

BED file contains the coordinates of the regions you'd like to keep from the downloaded VCF.

- Make sure the contig name in the BED file matches one of the contig names from the url list.
- If a contig name has more than one occurrence in the url list, its first value will be picked.

Last but not least, be sure that your BED is in the same coordinates as the downloaded VCF. For instance, the provided URLs (`gnomad_urls_download.csv`) point to hg38 vcf files.

