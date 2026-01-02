use std::fs::File;
use std::io::{self, BufReader, BufRead};
use std::collections::HashSet;


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
        fasta_records.push(
            Fasta { 
                header: std::mem::take(&mut cur_header), 
                sequence: std::mem::take(&mut cur_sequence) 
            }
        );
    }

    Ok(fasta_records)
}



pub struct Fastq {
    pub header: String,
    pub sequence: String,
    pub quality: String,
    pub compatable_transcripts: HashSet<String>,
}


pub fn read_fastq(path: &str) -> io::Result<Vec<Fastq>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    let mut fastq_records = Vec::new();
    let mut cur_header = String::new();
    let mut cur_sequence = String::new();
    let mut cur_quality = String::new();
    let mut state: usize = 0;

    for line in reader.lines() {
        let unwrapped_line = line?;
        let unwrapped_line = unwrapped_line.trim();

        match state {
            0 => {
                if let Some(rest) = unwrapped_line.strip_prefix('@') {
                    if !cur_header.is_empty() {
                        fastq_records.push(Fastq {
                            header: std::mem::take(&mut cur_header),
                            sequence: std::mem::take(&mut cur_sequence),
                            quality: std::mem::take(&mut cur_quality),
                            compatable_transcripts: HashSet::new(),
                        });
                    }
                    cur_header = rest.to_string();
                    state = 1;
                }
            },
            1 => {
                cur_sequence.push_str(unwrapped_line);
                state = 2;
            },
            2 => {
                if unwrapped_line.starts_with('+') {
                    state = 3;
                }
            },
            3 => {
                cur_quality.push_str(unwrapped_line);
                state = 0;
            },
            _ => {}
        }
    }

    if !cur_header.is_empty() {
        fastq_records.push(
            Fastq { 
                header: std::mem::take(&mut cur_header), 
                sequence: std::mem::take(&mut cur_sequence), 
                quality: std::mem::take(&mut cur_quality), 
                compatable_transcripts: HashSet::new() 
            }
        );
    }

    Ok(fastq_records)
}