use pyo3::prelude::*;
use pyo3::create_exception;
use pyo3::exceptions::PyException;

mod qualtrim;

create_exception!(qualtrim, HasNoQualities, PyException);


#[pyfunction(base=33)]
fn quality_trim_index(qualities: &str, cutoff_front: i32, cutoff_back: i32, base: u8) -> (usize, usize) {
    // TODO cutoff_front
    // TODO base
    (0, qualtrim::quality_trim_index(qualities.as_bytes(), cutoff_back))
}

#[pyfunction(base=33)]
fn nextseq_trim_index(sequence: &str, qualities: &str, cutoff: i32, base: u8) -> usize {
    // TODO base
    qualtrim::twocolor_trim_index(sequence.as_bytes(), qualities.as_bytes(), cutoff)
}

#[pyfunction(base=33)]
fn expected_errors(qualities: &str, base: u8) -> f64 {
    qualtrim::expected_errors(qualities.as_bytes(), base)
}

#[pymodule]
#[pyo3(name = "qualtrim")]
fn cutadapt_rust(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(quality_trim_index, m)?)?;
    m.add_function(wrap_pyfunction!(nextseq_trim_index, m)?)?;
    m.add_function(wrap_pyfunction!(expected_errors, m)?)?;
    m.add("HasNoQualities", py.get_type::<HasNoQualities>())?;
    Ok(())
}
