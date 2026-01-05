mod k_mer_utils;
mod file_parser_utils;
mod em_algorithm;
mod de_bruijn_graph;
mod file_writer_utils;

fn main() {
    let mut dbg = de_bruijn_graph::DeBruijnGraph::new(20);
    dbg.load_data("test_data/GJJC01.fasta", "test_data/SRR15662082_1.fastq");
    dbg.build_index_graph();
        
    let abundances = dbg.run_pseudoalignment(1000, 1e-3, 1e-6);

    file_writer_utils::abundance_to_csv(&abundances, "abundance_output.csv").expect("Failed to write abundance to CSV");
}
