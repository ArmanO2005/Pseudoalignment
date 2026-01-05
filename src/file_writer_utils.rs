use std::collections::HashMap;
use std::error::Error;
use csv::Writer;


pub fn abundance_to_csv(
    abundances: &HashMap<String, HashMap<String, f64>>,
    path: &str
) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(path)?;

    for (read_id, abundances) in abundances {
        for (transcript_id, abundance) in abundances {
            wtr.write_record(&[read_id, transcript_id, &abundance.to_string()])?;
        }
    }

    wtr.flush()?;
    Ok(())
} 