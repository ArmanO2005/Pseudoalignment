mod k_mer_utils;
mod file_parser_utils;
mod EM_algorithm;
mod de_bruijn_graph;

fn main() {
    let mut dbg = de_bruijn_graph::DeBruijnGraph::new(31);
    dbg.load_data("../test_data/GJJC01.fasta", "../test_data/SRR15662082_1.fastq");

    dbg!("data loaded");

    dbg.build_index_graph();

    dbg!("index built");
    
    let abundances = dbg.run_pseudoalignment(1000, 1e-6);
    
    for (read_id, read_abundance) in abundances {
        dbg!(read_id);
        dbg!(read_abundance);
    }
}
