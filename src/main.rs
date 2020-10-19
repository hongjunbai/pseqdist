mod pseqdist;

use std::env;
use std::fs::File;
use std::io::{Write, Error};
use std::path::Path;
use clap::{App, load_yaml};

fn main() -> Result<(), Error> {
    let yaml = load_yaml!("cli.yml");
    let args = App::from(yaml).get_matches();

    // Number of thread
    let max_cpus = num_cpus::get_physical() as i32;
    let mut nthread = args.value_of("nthread").unwrap()
                          .to_string().parse::<i32>().unwrap();
    if nthread <= 0 || nthread > max_cpus {
        nthread = max_cpus;
    } 
    // println!("nthread: {}", nthread);
    rayon::ThreadPoolBuilder::new().num_threads(nthread as usize).build_global().unwrap();

    // Input alignment
    // println!("Using input file: {}", args.value_of("INPUT").unwrap());
    let input = args.value_of("INPUT").unwrap().to_string();
    let (headers, seqs) = pseqdist::read_fasta(input);
    // for h in &headers { println!("{:?}", h); }
    let names: Vec<_> = headers.iter().map(|x| x.split(' ').next().unwrap()).collect();
    // for (n, s) in pseqdist::zip!(names, seqs) { println!(">{:?} | {:?}", n, s); }
    // Distance method
    let method = args.value_of("method").unwrap();
    let f_dist: fn(&str, &str) -> f32 = match method {
        "similarity" => pseqdist::similarity,
        "identity" => pseqdist::identity,
        "hamming" => hamming_f,
        _ => pseqdist::identity,
    };

    // Get the distance matrix
    // println!("Distance method: {}", args.value_of("method").unwrap());
    let pdist = pseqdist::pairwise_dist(&seqs, f_dist);
    let diag_fill = match method {
        "hamming" => 0.0,
        _ => 100.0,
    };
    let mat = pseqdist::to_mat(&pdist, diag_fill);

    // Output
    // println!("{:?}", mat);
    let outfile = args.value_of("outfile").unwrap().to_string();
    let path = Path::new(&outfile);
    let display = path.display();
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why),
        Ok(file) => file,
    };
    // Write
    // println!("Writing output file: {}", args.value_of("outfile").unwrap());
    let input_args: Vec<_> = env::args().collect();
    writeln!(file, "# {}", input_args.join(" "))?;
    writeln!(file, "#")?;
    writeln!(file, "seqs\t{}", names.join("\t"))?;
    for (name, row) in names.iter().zip(mat.genrows()) {
        let row_rec: Vec<_> = row.iter().map(|x| x.to_string()).collect();
        writeln!(file, "{}\t{}", name, row_rec.join("\t"))?;
    }

    Ok(())
}

// convert the return type as f32 to facilitate match
fn hamming_f(a: &str, b:&str) -> f32 {
    pseqdist::hamming(&a, &b) as f32
}
