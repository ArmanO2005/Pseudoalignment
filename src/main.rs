mod k_mer_utils;
mod file_parser_utils;
mod EM_algorithm;
mod de_bruijn_graph;

fn main() {
    let mut dbg = de_bruijn_graph::DeBruijnGraph::new(31);
    dbg.load_data("../data/transcripts.fasta", "../data/reads.fastq");
    dbg.build_index_graph();
    
    let abundances = dbg.run_pseudoalignment(1000, 1e-6);
    
    for (transcript_id, abundance_map) in abundances {
        println!("Transcript: {}", transcript_id);
        for (method, abundance) in abundance_map {
            println!("  {}: {:.6}", method, abundance);
        }
    }
}
