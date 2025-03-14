use std::ffi::{c_char, c_float, c_int, CStr, CString};

pub trait TSubread {
    fn get_seq(&self) -> &str;
    fn get_cx(&self) -> u8;
}

#[repr(C)]
pub struct Subread {
    seq: *const c_char,
    flags: c_int,
}

impl Subread {
    pub fn new(seq: *const c_char, flags: c_int) -> Self {
        Self { seq, flags }
    }
}

#[repr(C)]
pub struct Setting {
    min_identity: c_float,
    match_score: c_int,
    mismatch_score: c_int,
    insertion_score: c_int,
    deletion_score: c_int,
}

impl Default for Setting {
    fn default() -> Self {
        Self {
            min_identity: 0.82,
            match_score: 3,
            mismatch_score: -5,
            insertion_score: -2,
            deletion_score: -2,
        }
    }
}

#[repr(C)]
pub struct Result {
    seq: *mut c_char,
    n_passes: usize,
}

extern "C" {
    fn PoaDraftGen(reads: *const Subread, num_sbr: usize, setting: *const Setting) -> Result;
    fn FreePoaResult(result: Result);
}

pub fn poa_consensus<T: TSubread>(read_infos: &Vec<T>, setting: &Setting) -> (String, usize) {
    let read_infos = read_infos
        .iter()
        .map(|read_info| {
            (
                CString::new(read_info.get_seq()).unwrap(),
                read_info.get_cx() as c_int,
            )
        })
        .collect::<Vec<_>>();
    let reads = read_infos
        .iter()
        .map(|(seq, flag)| Subread::new(seq.as_ptr(), *flag as c_int))
        .collect::<Vec<_>>();

    unsafe {
        let poa_res = PoaDraftGen(reads.as_ptr(), reads.len(), setting);
        // let seq = CString::from_raw(poa_res.seq).to_str().unwrap().to_string();
        let seq = CStr::from_ptr(poa_res.seq).to_str().unwrap().to_string();
        let npasses = poa_res.n_passes;
        FreePoaResult(poa_res);
        (seq, npasses)
    }
}
