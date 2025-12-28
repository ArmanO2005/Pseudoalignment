mod k_mer_utils;
mod file_parser_utils;
use std::collections::{HashMap, HashSet};


struct DeBruijnGraph {
    pub k: usize,
    pub adj_list: HashMap<String, Vec<String>>,
    pub edge_index: HashMap<(String, String), HashSet<String>>,
}


impl DeBruijnGraph {
    pub fn new(k: usize) -> Self {
        assert!(k >= 2, "k must be >= 2");
        DeBruijnGraph {
            k,
            adj_list: HashMap::new(),
            edge_index: HashMap::new(),
        }
    }


    pub fn build_index_graph(&mut self, fasta_path: &str) {
        let fasta_records = file_parser_utils::read_fasta(fasta_path).expect("Failed to read fasta");
        let k_mer_index = k_mer_utils::build_k_mer_index(&fasta_records, self.k);

        for (transcript_id, k_mers) in k_mer_index.index {
            for i in 0..k_mers.len() {
                let k_mer = &k_mers[i];
                let k_prefix = &k_mer[..self.k - 1];
                let k_suffix = &k_mer[1..];

                let key = (k_prefix.to_string(), k_suffix.to_string());
                self.edge_index.entry(key.clone()).or_default().insert(transcript_id.clone());
                self.adj_list.entry(key.0.clone()).or_default().push(key.1.clone());
                self.adj_list.entry(key.1.clone()).or_default();
                

                self.adj_list.get_mut(k_prefix).unwrap().push(k_suffix.to_string());

            }
        }
    }


    fn pseudoalign(&self, read_sequence: &str) -> (HashSet<String>, Vec<(String, String)>) {
        let read_k_mers = k_mer_utils::create_k_mers(read_sequence, self.k);
        let mut cur: Option<HashSet<String>> = None;
        let mut read_edges = Vec::new();

        for k_mer in read_k_mers {
            let k_prefix = &k_mer[..self.k - 1];
            let k_suffix = &k_mer[1..];

            let key = (k_prefix.to_string(), k_suffix.to_string());

            read_edges.push(key.clone());

            let Some(mer) = self.edge_index.get(&key) else {
                continue;
            };

            match &mut cur {
                None => {
                    cur = Some(mer.iter().cloned().collect());
                }
                Some(cur_set) => {
                    cur_set.retain(|t| mer.contains(t));
                    if cur_set.is_empty() {
                        break;
                    }
                }
            }
        }

        (cur.unwrap_or_default(), read_edges)
    }

}
