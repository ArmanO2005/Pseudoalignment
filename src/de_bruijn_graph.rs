use std::collections::{HashMap, HashSet, BTreeMap};
use crate::file_parser_utils::{Fasta, Fastq};
use crate::k_mer_utils;
use crate::file_parser_utils;
use crate::em_algorithm;

pub struct DeBruijnGraph {
    pub k: usize,
    pub adj_list: HashMap<String, Vec<String>>,
    pub edge_index: HashMap<(String, String), HashSet<String>>,
    pub eff_lengths: HashMap<String, usize>,
    pub transcripts: Vec<Fasta>,
    pub reads: Vec<Fastq>,
}


impl DeBruijnGraph {
    pub fn new(k: usize) -> Self {
        assert!(k >= 2, "k must be >= 2");
        DeBruijnGraph {
            k,
            adj_list: HashMap::new(),
            edge_index: HashMap::new(),
            eff_lengths: HashMap::new(),
            transcripts: Vec::new(),
            reads: Vec::new(),
        }
    }


    pub fn load_data(&mut self, fasta_path: &str, fastq_path: &str) {
        self.transcripts = file_parser_utils::read_fasta(fasta_path).expect("Failed to read fasta");
        self.reads = file_parser_utils::read_fastq(fastq_path).expect("Failed to read fastq");
        // self.reads.truncate(30000);
    }


    pub fn build_index_graph(&mut self) {
        let k_mer_index = k_mer_utils::build_k_mer_index(&self.transcripts, self.k);

        for (transcript_id, k_mers) in k_mer_index.index {
            for i in 0..k_mers.len() {
                let k_mer = &k_mers[i];
                let k_prefix = &k_mer[..self.k - 1];
                let k_suffix = &k_mer[1..];

                let key = (k_prefix.to_string(), k_suffix.to_string());
                self.edge_index.entry(key.clone()).or_default().insert(transcript_id.clone());
                self.adj_list.entry(key.0.clone()).or_default().push(key.1.clone());
                self.adj_list.entry(key.1.clone()).or_default();
            }
        }

        for fasta in &self.transcripts {
            let eff_length = fasta.sequence.len().saturating_sub(self.k - 1);
            self.eff_lengths.insert(fasta.header.clone(), eff_length);
        }
    }


    fn get_compatability(&self, read_sequence: &str) -> (Vec<String>, Vec<(String, String)>) {
        let read_k_mers = k_mer_utils::create_k_mers(read_sequence, self.k);
        let mut cur: Option<HashSet<String>> = None;
        let mut read_edges = Vec::new();

        for k_mer in read_k_mers {
            let k_prefix = &k_mer[..self.k - 1];
            let k_suffix = &k_mer[1..];

            let key = (k_prefix.to_string(), k_suffix.to_string());

            read_edges.push(key.clone());

            let Some(transcript_list) = self.edge_index.get(&key) else {
                return (Vec::new(), read_edges);
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

        let result = cur.unwrap_or_default();
        (result.into_iter().collect::<Vec<_>>(), read_edges)
    }


    fn get_equivalence_classes(
        &mut self, 
    ) -> BTreeMap<Vec<String>, usize> {
        let mut eq_classes: BTreeMap<Vec<String>, usize> = BTreeMap::new();
        
        let mut compat_data = Vec::new();
        for record in self.reads.iter() {
            let (compatability_set, _) = self.get_compatability(&record.sequence);
            compat_data.push(compatability_set);
        }
        
        for (i, compatability_set) in compat_data.into_iter().enumerate() {
            let mut sorted_set = compatability_set.clone();
            sorted_set.sort();
            *eq_classes.entry(sorted_set.clone()).or_insert(0) += 1;
            self.reads[i].compatable_transcripts = compatability_set.into_iter().collect();
        }
        eq_classes
    }


    fn log_likelihood_fxn(&self, 
        alpha_param: HashMap<String, f64>, 
        eq_classes: &BTreeMap<Vec<String>, usize>, 
    ) -> f64 {

        assert!(alpha_param.len() == self.transcripts.len(), "Alpha parameter length must match number of transcripts");
        
        let mut log_likelihood = 0.0;

        for (transcript_ids, count) in eq_classes {
            let mut base: f64 = 0.0;
            for transcript_id in transcript_ids {
                let transcript_length = self.eff_lengths[transcript_id] as f64;
                base += alpha_param[transcript_id] / transcript_length;
            }
            log_likelihood += *count as f64 * base.ln();
        }
        log_likelihood
    }


    fn approx_MLE(&mut self, 
        max_iterations: usize,
        negligable_error: f64
    ) -> HashMap<String, f64> {
        let eq_classes = self.get_equivalence_classes();

        let log_likelihood_fn = |alpha: HashMap<String, f64>, 
                             eq_classes: &BTreeMap<Vec<String>, usize>| {
            self.log_likelihood_fxn(alpha, eq_classes)
        };

        em_algorithm::expectation_maximization_algorithm(
            &log_likelihood_fn,
            &self.transcripts,
            &eq_classes,
            self.k,
            max_iterations,
            negligable_error,
        )
    }


    fn get_abundances(&self, approx_MLE: HashMap<String, f64>) -> HashMap<String, HashMap<String, f64>> {

        let mut abundances: HashMap<String, HashMap<String, f64>> = HashMap::new();
        
        for read in self.reads.iter() {
            let transcripts = read.compatable_transcripts.iter().collect::<Vec<_>>();
            let mut read_abundances: HashMap<String, f64> = HashMap::new();

            let denominator = transcripts.iter().map(|x| approx_MLE[*x]).sum::<f64>();

            for transcript_id in transcripts {
                let prob = approx_MLE[transcript_id] / denominator;
                read_abundances.insert(transcript_id.clone(), prob);
            }

            abundances.insert(read.header.clone(), read_abundances);
        } 

        abundances
    }


    fn prune_transcripts(&self, 
        abundances: &HashMap<String, HashMap<String, f64>>,
        tolerance: f64
    ) -> HashMap<String, HashMap<String, f64>> {
        let mut pruned_abundances: HashMap<String, HashMap<String, f64>> = HashMap::new();

        for (read_id, read_abundances) in abundances {
            for (transcript_id, abundance) in read_abundances {
                if *abundance >= tolerance {
                    pruned_abundances.entry(read_id.clone()).or_default().insert(transcript_id.clone(), *abundance);
                }
            }
        }

        pruned_abundances
    }

    pub fn run_pseudoalignment(&mut self, 
        max_iterations: usize,
        negligable_error: f64,
        tolerance: f64
    ) -> HashMap<String, HashMap<String, f64>> {
        let mle = self.approx_MLE(max_iterations, negligable_error);
        let abundances = self.get_abundances(mle);
        let pruned_abundances = self.prune_transcripts(&abundances, tolerance);

        pruned_abundances
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

}
