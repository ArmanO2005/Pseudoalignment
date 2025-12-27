use std::fs::File;
use std::io::{self, BufReader, BufRead};


pub struct Fasta {
    pub header: String,
    pub sequence: String,
}


pub fn read_fasta(path: &str) -> io::Result<Vec<Fasta>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    let mut fasta_records = Vec::new();
    let mut cur_header = String::new();
    let mut cur_sequence = String::new();


    for line in reader.lines() {
        let unwrapped_line = line?;
        let unwrapped_line = unwrapped_line.trim();

        if let Some(rest) = unwrapped_line.strip_prefix('>') {
            if !cur_header.is_empty() {
                fasta_records.push(Fasta {
                    header: std::mem::take(&mut cur_header),
                    sequence: std::mem::take(&mut cur_sequence),
                });
            }
            cur_header = rest.to_string();
        } else {
            cur_sequence.push_str(unwrapped_line);
        }
    }

    if !cur_header.is_empty() {
        fasta_records.push(Fasta { header: std::mem::take(&mut cur_header), sequence: std::mem::take(&mut cur_sequence) });
    }

    Ok(fasta_records)
}



pub struct Fastq {
    pub header: String,
    pub sequence: String,
    pub quality: String,
}


pub fn read_fastq(path: &str) -> io::Result<Vec<Fastq>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    let mut fastq_records = Vec::new();
    let mut cur_header = String::new();
    let mut cur_sequence = String::new();
    let mut cur_quality = String::new();
    let mut on_seq: bool = true;

    for line in reader.lines() {
        let unwrapped_line = line?;
        let unwrapped_line = unwrapped_line.trim();

        if let Some(rest) = unwrapped_line.strip_prefix('@') {
            if !cur_header.is_empty() {
                fastq_records.push(Fastq {
                    header: std::mem::take(&mut cur_header),
                    sequence: std::mem::take(&mut cur_sequence),
                    quality: std::mem::take(&mut cur_quality),
                });
            }
            cur_header = rest.to_string();
        } else if unwrapped_line.starts_with('+') {
            on_seq = false;
            continue;       
        } else if on_seq{
            cur_sequence.push_str(unwrapped_line)
        } else {
            cur_quality.push_str(unwrapped_line);
            on_seq = true;
        }
    }

    if !cur_header.is_empty() {
        fastq_records.push(Fastq { header: std::mem::take(&mut cur_header), sequence: std::mem::take(&mut cur_sequence), quality: std::mem::take(&mut cur_quality) });
    }

    Ok(fastq_records)
}