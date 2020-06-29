use std::fs::File;
use std::io::{self, BufRead};
use std::collections::HashSet;

//extern crate itertools;
use itertools::Itertools;
use rayon::prelude::*;
use ndarray::prelude::*;
//use std::thread;

#[macro_export]
macro_rules! zip { ($a:ident, $b:ident) => { $a.iter().zip($b.iter()) } }
macro_rules! zip_str { ($a:ident, $b:ident) => { $a.chars().zip($b.chars()) } }

pub fn read_lines(filename: std::string::String) -> io::Lines<io::BufReader<File>> {
    // Returns an Iterator to the Reader of the lines of the file.
    let file = match File::open(&filename) {
        Err(why) => panic!("couldn't open {}: {}", &filename, why),
        Ok(file) => file,
    };
    io::BufReader::new(file).lines()
}

pub fn read_fasta (infile: std::string::String) -> (Vec<String>, Vec<String>) {
    let mut headers: Vec<String> = vec![];
    let mut sequences: Vec<String> = vec![];
    for line in read_lines(infile) {
        let cleaned = line.unwrap().trim().to_string(); //replace(" ", "");
        // println!("{:?}", cleaned);
        if cleaned.starts_with('>') {
            headers.push(cleaned[1..].trim().to_string());
            sequences.push("".to_string())
        } else {
            sequences.last_mut().unwrap().push_str(&cleaned.replace(" ", ""))
        }
    }
    (headers, sequences)
}

pub fn hamming(a: &str, b: &str) -> usize {
    zip_str!(a, b).filter(|&(char_a, char_b)| char_a != char_b).count()
}

pub fn identity(a: &str, b: &str) -> f32 {
    // Gap-compressed identity: consecutive gaps count as one difference
    let mut matched_len = 0.0;
    let mut total_len = 0.0;
    let mut is_indel = false;
    for (ai, bi) in zip_str!(a, b) {
        if ai == '-' && bi == '-' { // empty column
            continue;
        } else if ai == '-' || bi == '-' { // gap 
            // consecutive gaps count as one difference
            if ! is_indel {
                is_indel = true;
                total_len +=1.0;
            }
        } else if ai == bi {  // matched
            matched_len += 1.0;
            total_len += 1.0;
            is_indel = false;
        } else { // mismatched, non-gap
            total_len += 1.0;
            is_indel = false;
        }
    }
    matched_len * 100.0 / total_len
}

pub fn similarity(a: &str, b: &str) -> f32 {
    let similar: HashSet<(_, _)> = "RK KR DE ED ND DN QE EQ QN NQ ST TS SA AS VI IV IL LI LM ML FY YF"
                                   .split(' ')
                                   .map(|p| {let v: Vec<_> = p.chars().collect(); (v[0], v[1])})
                                   .collect();
    // Gap-compressed similarity: consecutive gaps count as one difference
    let mut matched_len = 0.0;
    let mut total_len = 0.0;
    let mut is_indel = false;
    for (ai, bi) in zip_str!(a, b) {
        if ai == '-' && bi == '-' { // empty column
            continue;
        } else if ai == '-' || bi == '-' { // gap 
            // consecutive gaps count as one difference
            if ! is_indel {
                is_indel = true;
                total_len +=1.0;
            }
        } else if ai == bi || similar.contains(&(ai, bi)) {  // matched
            matched_len += 1.0;
            total_len += 1.0;
            is_indel = false;
        } else { // mismatched, non-gap
            total_len += 1.0;
            is_indel = false;
        }
    }
    matched_len * 100.0 / total_len
}

pub fn pairwise_dist<T: Send> (seqs: &[String], f: fn(&str, &str) -> T) -> Vec<T> {
    // parallel using rayon (.par_iter())
    let pidxs: Vec<_> = (0..seqs.len()).tuple_combinations::<(_,_)>().collect();
    pidxs.par_iter().map(|(i, j)| f(&seqs[*i as usize], &seqs[*j as usize])).collect()
}

pub fn to_mat<T: Clone+Copy> (pdist: &[T], diag_fill: T) -> Array2<T> {
    let t = (1.0 + ((1+8*pdist.len()) as f64).sqrt())/2.0;
    assert!((t-t.round()).abs() < 1E-6);  // Fail if pdist is of wrong size (didn't correspond to a matrix)
    let n = t as usize;
    //
    let mut result = Array::from_elem((n, n), diag_fill);
    let mut ii = 0;
    for i in 0..n {
        for j in i+1..n {
            result[[i,j]] = pdist[ii];
            result[[j,i]] = pdist[ii];
            ii += 1;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! vec_of_strings { ($($x:expr),*) => (vec![$($x.to_string()),*]); }

    #[test]
    fn test_read_fasta() {
        let expected_names = vec_of_strings!["seq1", "seq2", "seq3"];
        let expected_seqs = vec_of_strings!["AD-DSE-FPHTTSTST", "AD-DSE-FPHF-SS-H", "ADLHDH-SGHA-SS--"];
        let (names, seqs) = read_fasta("tests/data/toy.fa".to_string());
        let n_matching_names = zip!(expected_names, names).filter(|&(a, b)| a == b).count();
        let n_matching_seqs = zip!(expected_seqs, seqs).filter(|&(a, b)| a == b).count();
        assert!(n_matching_names == expected_names.len());
        assert!(n_matching_seqs == expected_seqs.len());
    }

    #[test]
    fn test_identity() {
        let seq1 = "AD---E-FPHTTSTST".to_string();
        let seq2 = "AD-DSE-FPHF--SSH".to_string();
        let expected_identity = 7.0 * 100.0/12.0;
        let identity = identity(&seq1, &seq2);
        println!("{:?}", identity);
        assert!(expected_identity == identity);
    }

    #[test]
    fn test_similarity() {
        let seq1 = "AD---E-FPHVTSTSY".to_string();
        let seq2 = "AD-DSE-FPHI--SSF".to_string();
        let expected_similarity = 10.0 * 100.0/12.0;
        let sim = similarity(&seq1, &seq2);
        println!("{:?}", sim);
        assert!(expected_similarity == sim);
    }

    #[test]
    fn test_hamming() {
        let seq1 = "AD-DSE-FPHTTSTST".to_string();
        let seq2 = "AD-DSE-FPHF-SS-H".to_string();
        let expected_hamming = 5;
        let hamming = hamming(&seq1, &seq2);
        assert!(expected_hamming == hamming);
    }

    #[test]
    fn test_pairwise_dist() {
        let seqs = vec_of_strings!["AD-DSE-FPHTTSTST", "AD-DSE-FPHF-SS-H", "ADLHDH-SGHA-SS--"];
        let expected_pdist = vec![5, 11, 8];
        let pdist = pairwise_dist(&seqs, hamming);
        let n_matches = zip!(expected_pdist, pdist).filter(|&(a, b)| a == b).count();
        assert!(n_matches == 3);
    }

    #[test]
    fn test_to_mat() {
        let pdist = vec![5, 11, 8];
        let mat = to_mat(&pdist, 0);
        let expected_mat = array![[ 0,  5, 11],
                                  [ 5,  0,  8],
                                  [11,  8,  0]];
        let n_matches = zip!(expected_mat, mat).filter(|&(a, b)| a == b).count();
        assert!(n_matches == 9);
    }
}
