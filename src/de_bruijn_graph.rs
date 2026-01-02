mod k_mer_utils;
mod file_parser_utils;
mod EM_algorithm;
use std::collections::{HashMap, HashSet};
use crate::file_parser_utils::{Fasta, Fastq};

fn main() {
    let fasta_records = file_parser_utils::read_fasta("../test_data/GJJC01.fasta").expect("Failed to read fasta");
    let fastq_records = file_parser_utils::read_fastq("../test_data/SRR15662082_1.fastq").expect("Failed to read fastq");
    let k = 7;
    let mut de_bruijn_graph = DeBruijnGraph::new(k);
    de_bruijn_graph.build_index_graph(&fasta_records);
    let scores = de_bruijn_graph.get_MLEs(
        fastq_records,
        fasta_records,
        100,
        1e-6,
    );
    dbg!(scores);
}


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


    pub fn build_index_graph(&mut self, fasta_records: &Vec<Fasta>) {
        let k_mer_index = k_mer_utils::build_k_mer_index(&fasta_records, self.k);

        for (transcript_id, k_mers) in k_mer_index.index {
            for i in 0..k_mers.len() {
                let k_mer = &k_mers[i];
                let k_prefix = &k_mer[..self.k - 1];
                let k_suffix = &k_mer[1..];

                let key = (k_prefix.to_string(), k_suffix.to_string());
                self.edge_index.entry(key.clone()).or_default().push(transcript_id.clone());
                self.adj_list.entry(key.0.clone()).or_default().push(key.1.clone());
                self.adj_list.entry(key.1.clone()).or_default();
            }
        }
    }


    fn get_compatability(&self, read_sequence: &str) -> (HashSet<String>, Vec<(String, String)>) {
        let read_k_mers = k_mer_utils::create_k_mers(read_sequence, self.k);
        let mut cur: Option<HashSet<String>> = None;
        let mut read_edges = Vec::new();

        for k_mer in read_k_mers {
            let k_prefix = &k_mer[..self.k - 1];
            let k_suffix = &k_mer[1..];

            let key = (k_prefix.to_string(), k_suffix.to_string());

            read_edges.push(key.clone());

            let Some(transcript_list) = self.edge_index.get(&key) else {
                return (HashSet::new(), read_edges);
            };

            match &mut cur {
                None => {
                    cur = Some(transcript_list.iter().cloned().collect());
                }
                Some(cur_set) => {
                    cur_set.retain(|t| transcript_list.contains(t));
                    if cur_set.is_empty() {
                        break;
                    }
                }
            }
        }

        (cur.unwrap_or_default(), read_edges)
    }


    fn get_equivalence_classes(&self, fastq_records: Vec<Fastq>) -> (HashMap<HashSet<String>, usize>, HashMap<String, HashSet<String>>) {
        let mut eq_classes: HashMap<HashSet<String>, usize> = HashMap::new();
        let mut read_compatabilities: HashMap<String, HashSet<String>> = HashMap::new();

        for record in fastq_records {
            let (compatability_set, _) = self.get_compatability(&record.sequence);
            *eq_classes.entry(compatability_set).or_default() += 1;
            read_compatabilities.entry(record.clone()).or_default().push(compatability_set);
        }
        (eq_classes, read_compatabilities)
    }


    fn log_likelihood_fxn(&self, 
        alpha_param: HashMap<String, f64>, 
        eq_classes: &HashMap<HashSet<String>, usize>, 
        fasta_records: Vec<Fasta>,
        eff_lengths: HashMap<String, f64>
    ) -> f64 {

        assert!(alpha_param.len() == fasta_records.len(), "Alpha parameter length must match number of transcripts");
        
        let mut log_likelihood = 0.0;

        for e in eq_classes {
            let mut base = 0.0;
            for transcript_id in &e.0 {
                let transcript_length = eff_lengths[transcript_id];
                base += alpha_param[transcript_id] / transcript_length;
            }
            log_likelihood += e.1 as f64 * base.ln();
        }
        log_likelihood
    }


    fn get_abundances(&self, 
        approx_MLE: HashMap<String, f64>, 
        read_compatabilities: HashMap<String, HashSet<String>>
    ) -> HashMap<String, HashMap<String, f64>> {
        let mut abundances: HashMap<String, HashMap<String, f64>> = HashMap::new();
        
        for read in read_compatabilities.index {
            let transcripts = read_compatabilities[read];
            let mut read_abundances: HashMap<String, f64> = HashMap::new();

            let denominator = transcripts.iter().map(|x| approx_MLE[x]).sum();

            for transcript_id in transcripts {
                let prob = approx_MLE[transcript_id] / denominator;

                read_abundances.entry(transcript_id.clone()).or_default() = prob;
            }

            abundances.entry(read.clone()).or_default() = read_abundances;
        } 

        abundances
    }


    fn approx_MLE(&self, 
        fastq_records: Vec<Fastq>, 
        fasta_records: Vec<Fasta>,
        max_iterations: usize,
        tolerance: f64
    ) -> HashMap<String, f64> {
        let eq_classes = self.get_equivalence_classes(fastq_records);

        let log_likelihood_fn = |alpha: HashMap<String, f64>, 
                             eq_classes: &HashMap<HashSet<String>, usize>,
                             fasta_records: &Vec<Fasta>,
                             eff_lengths: &HashMap<String, f64>| {
            self.log_likelihood_fxn(alpha, eq_classes, fasta_records.clone(), eff_lengths.clone())
        };

        EM_algorithm::expectation_maximization_algorithm(
            &log_likelihood_fn,
            &eq_classes,
            &fasta_records,
            self.k,  // Pass k!
            max_iterations,
            tolerance,
        )
    }






    // fn walk_score(&self, 
    //     read_edges: &[(String, String)], 
    //     eligable_transcripts: &HashSet<String>
    // ) -> HashMap<String, f64> {

    //     let mut transcript_scores = HashMap::new();
    //     for transcript_id in eligable_transcripts {
    //         let mut successes = Vec::new();
    //         for i in 0..read_edges.len() - 1 {
    //             let (_, suffix) = &read_edges[i];
    //             let (_, next_suffix) = &read_edges[i + 1];

    //             let key = (suffix.clone(), next_suffix.clone());

    //             if let Some(neighbors) = self.adj_list.get(suffix) {
    //                 if !neighbors.contains(next_suffix) {
    //                     successes.push(false);
    //                 } else if let Some(transcripts) = self.edge_index.get(&key) {
    //                     successes.push(transcripts.contains(transcript_id));
    //                 } else {
    //                     successes.push(false);
    //                 }
    //             } else {
    //                 successes.push(false);
    //             }
    //         }
    //         let score = successes.iter().filter(|&&x| x).count() as f64 / successes.len().max(1) as f64;            
    //         transcript_scores.insert(transcript_id.clone(), score);
    //     }
    //     transcript_scores
    // }


    // fn align_read(&self, read_sequence: &str) -> HashMap<String, f64> {
    //     let (eligable_transcripts, read_edges) = self.pseudoalign(read_sequence);
    //     self.walk_score(&read_edges, &eligable_transcripts)
    // }

    // pub fn align_fastq(&self, fastq_path: &str) -> Vec<HashMap<String, f64>> {
    //     let fastq_records = file_parser_utils::read_fastq(fastq_path).expect("Failed to read fastq");
    //     let mut results = Vec::new();
    //     let mut i = 0;
    //     for record in fastq_records {
    //         if i < 1000 {
    //             let alignment = self.align_read(&record.sequence);
    //             results.push(alignment);
    //         }
    //         i += 1;
    //     }
    //     results
    // }
}
