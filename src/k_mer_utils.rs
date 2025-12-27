use std::collections::HashMap;


pub fn create_k_mers(sequence: &str, k: usize) -> Vec<String> {
    let mut k_mers = Vec::new();
    let num_k_mers = sequence.len().saturating_sub(k - 1);

    for i in 0..num_k_mers {
        let k_mer = &sequence[i..i + k];
        k_mers.push(k_mer.to_string());
    }   

    k_mers
}


struct Fasta {
    pub header: String,
    pub sequence: String,
}


struct transcript_k_mer {
    pub transcript_id: String,
    pub k_mer: Vec<String>,
}

struct KMerIndex {
    pub index: HashMap<String, Vec<String>>,
}


pub fn build_k_mer_index(sequences: &Vec<Fasta>, k: usize) -> KMerIndex {
    let mut k_mer_index = KMerIndex {
        index: HashMap::new(),
    };

    for fasta in sequences {
        let k_mers = create_k_mers(&fasta.sequence, k);
        k_mer_index.index.insert(fasta.header.clone(), k_mers);
    }

    k_mer_index
}