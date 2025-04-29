# Why smart_gnomad_downloader?

This simple rust project, with no dependencies, allows one to download big VCF files such as gnomad, while only saving to disk what intersects a given BED file.

I know you can simply `wget` or `curl` the file and pipe it to `bedtools`, but this can be tedious, especially if you need to ensure the checksum matches.

# How to install

