const QUALITY_BASE: u8 = 33;

// Find the positions at which to trim low-quality ends from a nucleotide sequence.
// Return tuple (start, stop) that indicates the good-quality segment.
//
// Qualities are assumed to be ASCII-encoded as chr(qual + base).
//
// The algorithm is the same as the one used by BWA within the function
// 'bwa_trim_read':
// - Subtract the cutoff value from all qualities.
// - Compute partial sums from all indices to the end of the sequence.
// - Trim sequence at the index at which the sum is minimal.
pub fn quality_trim_index(qualities: &[u8], cutoff: i32) -> usize {
    let mut s = 0;
    let mut max_qual = 0;
    let mut end = qualities.len();

    for i in (0..qualities.len()).rev() {
        // println!("in qualtrimindex: i={}, s={}, cutoff={}, q[i]={}", i, s, cutoff, qualities[i]);
        s += cutoff - (qualities[i] - QUALITY_BASE) as i32;
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            end = i;
        }
    }
    end
}

/// Quality trimming that works with Illumina two-color (2-channel) chemistry
/// as introduced in NextSeq where "no color" (dark cycle) encodes a G base.
/// The issue is that an apparent run of high-quality Gs towards the end of a read is
/// likely to indicate that actually the end of the fragment has been reached.
///
/// This function uses a modified trimming algorithm that forces quality values
/// of all G bases to cutoff - 1.
pub fn twocolor_trim_index(sequence: &[u8], qualities: &[u8], cutoff: i32) -> usize {
    let mut s = 0;
    let mut max_qual = 0;
    let mut end = qualities.len();
    assert_eq!(sequence.len(), qualities.len());
    for i in (0..end).rev() {
        let q = if sequence[i] == b'G' {
            cutoff - 1
        } else {
            (qualities[i] - QUALITY_BASE) as i32
        };
        s += cutoff - q as i32;
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            end = i
        }
    }
    end
}

// Return the number of expected errors (as double) from a read’s
// qualities.
//
// This uses a formula from Edgar et al. (2015),
// see Section 2.2 in <https://academic.oup.com/bioinformatics/article/31/21/3476/194979>.
//
pub fn expected_errors(qualities: &[u8], base: u8) -> f64 {
    let mut e = 0.0;
    for i in 0..qualities.len() {
        let q = qualities[i] - base;
        e += 10f64.powf(-(q as f64) / 10.0);
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;
    /*
        fn test test_nextseq_trim() {
            assert_eq!(nextseq_trim_index(SequenceRecord
            s = Sequence("n", "", "")
            assert nextseq_trim_index(s, cutoff=22) == 0
            s = Sequence(
                "n",
                "TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG",
                "AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA",
            )
            assert nextseq_trim_index(s, cutoff=22) == 33
        }
    */

    #[test]
    fn test_expected_errors() {
        assert_eq!(expected_errors(b""), 0.0);
        assert_eq!(expected_errors(&vec![10u8 + 33u8]), 0.1);
        assert_eq!(expected_errors(&vec![20u8 + 33u8]), 0.01);
        assert_eq!(
            expected_errors(&vec![10u8 + 33u8, 20u8 + 33u8, 30u8 + 33u8]),
            0.111
        );
    }

    #[test]
    fn test_twocolor_quality_trim_index() {
        assert_eq!(twocolor_trim_index(b"", b"", 22), 0);
        let sequence = b"TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG";
        let qualities = b"AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA";
        assert_eq!(twocolor_trim_index(sequence, qualities, 22), 33);
    }
}
