use std::collections::{HashMap, HashSet};
use crate::file_parser_utils::Fasta;



pub fn expectation_maximization_algorithm(
    log_likelihood_fxn: &dyn Fn(
        HashMap<String, f64>, 
        &HashMap<HashSet<String>, usize>, 
        &Vec<Fasta>,
        &HashMap<String, f64>
    ) -> f64,
    eq_classes: &HashMap<HashSet<String>, usize>,
    fasta_records: &Vec<Fasta>,
    k: usize,
    max_iterations: usize,
    tolerance: f64,
) -> HashMap<String, f64> {

    let mut alpha: HashMap<String, f64> = fasta_records.iter().map(|f| (f.header.clone(), 1.0)).collect();
    
    let eff_lengths: HashMap<String, f64> = fasta_records.iter().map(|f| (f.header.clone(), (f.sequence.len() - k + 1) as f64)).collect();
    
    let mut prev_log_likelihood = f64::NEG_INFINITY;
    
    for iteration in 0..max_iterations {

        let mut expected_counts: HashMap<String, f64> = fasta_records.iter().map(|f| (f.header.clone(), 0.0)).collect();
        
        for (transcript_set, read_count) in eq_classes {
            if transcript_set.is_empty() {
                continue;
            }
            
            let mut denominator = 0.0;
            for transcript_id in transcript_set {
                denominator += alpha[transcript_id] / eff_lengths[transcript_id];
            }
            
            if denominator == 0.0 {
                continue;
            }
            
            for transcript_id in transcript_set {
                let numerator = alpha[transcript_id] / eff_lengths[transcript_id];
                let probability = numerator / denominator;
                
                expected_counts.get_mut(transcript_id).map(|count| {*count += probability * (*read_count as f64)});
            }
        }
        

        alpha = expected_counts;
        
        let total: f64 = alpha.values().sum();
        if total > 0.0 {
            for value in alpha.values_mut() {
                *value /= total;
            }
        }
        
        let current_log_likelihood = log_liklihood_fxn(
            alpha.clone(), 
            eq_classes, 
            fasta_records,
            &eff_lengths
        );
        
        let likelihood_change = (current_log_likelihood - prev_log_likelihood).abs();
        
        
        if likelihood_change < tolerance {
            println!("Converged after {} iterations", iteration + 1);
            break;
        }
        
        prev_log_likelihood = current_log_likelihood;
    }
    
    alpha
}