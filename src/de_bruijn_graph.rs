mod k_mer_utils;
mod file_parser_utils;
use std::collections::{HashMap, HashSet};


struct DeBruijnGraph {
    pub k: usize,
    pub adj_list: HashMap<String, Vec<String>>,
    pub transcript_index: HashMap<String, HashSet<String>>,
}


impl DeBruijnGraph {
    fn new(k: usize) -> Self {
        DeBruijnGraph {
            k,
            adj_list: HashMap::new(),
            transcript_index: HashMap::new(),
        }
    }


    fn build_index_graph(&mut self, fasta_path: &str) {
        let fasta_records = file_parser_utils::read_fasta(fasta_path).expect("Failed to read fasta");
        let k_mer_index = k_mer_utils::build_k_mer_index(&fasta_records, self.k);

        for (transcript_id, k_mers) in k_mer_index.index {
            for i in 0..k_mers.len() {
                let k_mer = &k_mers[i];
                let k_prefix = &k_mer[..self.k - 1];
                let k_suffix = &k_mer[1..];

                if !self.adj_list.contains_key(k_prefix) {
                    self.adj_list.insert(k_prefix.to_string(), Vec::new());
                }

                if !self.adj_list.contains_key(k_suffix) {
                    self.adj_list.insert(k_suffix.to_string(), Vec::new());
                }

                if !self.transcript_index.contains_key(k_prefix) {
                    self.transcript_index.insert(k_prefix.to_string(), HashSet::new());
                }
                self.transcript_index.get_mut(k_prefix).unwrap().insert(transcript_id.clone());


                if !self.transcript_index.contains_key(k_suffix) {
                    self.transcript_index.insert(k_suffix.to_string(), HashSet::new());
                }
                self.transcript_index.get_mut(k_suffix).unwrap().insert(transcript_id.clone());


                self.adj_list.get_mut(k_prefix).unwrap().push(k_suffix.clone());

            }
        }
    }
}